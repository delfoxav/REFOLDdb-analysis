import pandas as pd
import textwrap
from bs4 import BeautifulSoup
import requests
from tqdm import tqdm
import time
import asyncio
import aiohttp
from datetime import datetime
import sqlite3
from urllib.parse import urljoin
from tqdm.asyncio import tqdm_asyncio
from aiohttp import ClientError, ClientConnectorError
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np

semaphore = asyncio.Semaphore(10)  # max 10 concurrent requests

MonashURLs = []
DOIs = []
session = requests.Session()

BaseURL = "https://pford.info/refolddb/"
search_all = "./list.cgi?param=*&lang=EN"
#search_all = "./list.cgi?kind=40&start=&end=null&lang=EN"
# first we get the col_names from a single record
r = session.get("https://pford.info/refolddatabase/refoldingrecord/5/")

soup = BeautifulSoup(r.text, "html.parser")
header = ""
col_names = []
data = []

for th in soup.find_all("th"):
    if th.get("class")[0] == "detail_header":
        header = th.find(string=True, recursive=False).strip()
    else:
        col_names.append(header + "." + th.text.replace(" ","_"))

col_names.append("DOI")
print("Gathered all the col names")
check_paper_row = ["Check the paper"] * (len(col_names) - 1)
no_data_row = ["NoData"]*len(col_names)

search_all_url = urljoin(BaseURL,search_all)

r = session.get(search_all_url)

soup = BeautifulSoup(r.text, "html.parser")
refs = soup.find_all("a",string="Detail")

links = [urljoin(BaseURL,link.get("href")) for link in refs]

print("Got all the records URLs, we have ", len(links), "records")

refoldDF = pd.read_table("REFOLD_161028.txt")

refoldDF["pubmed_id"] = refoldDF["pubmed_id"].astype("Int64").astype("str")
print("Loaded RefoldDF")

# Function to get true columns in a row
def true_columns(row):
    return ', '.join(row.index[row.astype(bool)])

async def FindMonashUrl(url, session):
    
    async with semaphore:
        try:
            async with session.get(url, timeout=aiohttp.ClientTimeout(total=10)) as response:
                if response.status == 429:
                    retry_after = int(response.headers.get("Retry-After", 5))
                    print(f"Rate limited, sleeping {retry_after} seconds...")
                    await asyncio.sleep(retry_after)
                response.raise_for_status()  # raise error on bad status
                html_content = await response.text()
        except (aiohttp.ClientError, asyncio.TimeoutError) as e:
            print(f"Error fetching {url}: {e}")
            return None
   
        soup = BeautifulSoup(html_content, "html.parser")

        # Find all the table rows
        rows = soup.find_all("tr")

        for row in rows:
            # checks if one of the column is "Link"
            if row.find("td", string="Link") != None:
                # get the url
                return urljoin(BaseURL,row.find("a").get("href"))
    
    return None

async def FindPubMedURL(url, session):

    async with semaphore:
        try:
            async with session.get(url, timeout=aiohttp.ClientTimeout(total=10)) as response:
                if response.status == 429:
                    retry_after = int(response.headers.get("Retry-After", 5))
                    print(f"Rate limited, sleeping {retry_after} seconds...")
                    await asyncio.sleep(retry_after)
                response.raise_for_status()  # raise error on bad status
                html_content = await response.text()
        except (aiohttp.ClientError, asyncio.TimeoutError) as e:
            print(f"Error fetching {url}: {e}")
            return None

        soup = BeautifulSoup(html_content, "html.parser")

        # Find all the table rows
        rows = soup.find_all("tr")

        for row in rows:
            # checks if one of the column is "PubMed ID"
            if row.find("td", string="PubMed ID") != None:
                # get the url
                return row.find("a").get("href")
        
    return None

async def FindDOI(url, session):

    async with semaphore:
        try:
            async with session.get(url, timeout=aiohttp.ClientTimeout(total=10)) as response:
                if response.status == 429:
                    retry_after = int(response.headers.get("Retry-After", 5))
                    print(f"Rate limited, sleeping {retry_after} seconds...")
                    await asyncio.sleep(retry_after)
                response.raise_for_status()  # raise error on bad status
                html_content = await response.text()
        except (aiohttp.ClientError, asyncio.TimeoutError) as e:
            print(f"Error fetching {url}: {e}")
            return None

    soup = BeautifulSoup(html_content, "html.parser")

    return soup.find("meta", attrs={"name": "citation_doi"}).get("content")

async def gatherrecord(url, session):
    """
    Function to gather all the value of a single record
    """
    async with semaphore:
        try:
            async with session.get(url, timeout=aiohttp.ClientTimeout(total=10)) as response:
                if response.status == 429:
                    retry_after = int(response.headers.get("Retry-After", 5))
                    print(f"Rate limited, sleeping {retry_after} seconds...")
                    await asyncio.sleep(retry_after)
                response.raise_for_status()  # raise error on bad status
                html_content = await response.text()
        except (aiohttp.ClientError, asyncio.TimeoutError) as e:
            print(f"Error fetching {url}: {e}")
            return None
    
        soup = BeautifulSoup(html_content, "html.parser")

        data = []
        for row in soup.find_all('tr'):
                for cell in row.find_all("td"):
                    data.append(cell.text)
        
        return data

async def process_url(url, session):
    row_record = None
    monash_url = await FindMonashUrl(url, session)
    if monash_url is not None:
        row_record = await gatherrecord(monash_url, session)
        row_record.append("Not needed")
        return row_record
    elif monash_url is None:
        row_record = await extract_from_db(url, session)
    pubmed_url = await FindPubMedURL(url, session)
    if pubmed_url is not None:
        if not row_record:
            row_record = check_paper_row.copy()
        DOI = await FindDOI(pubmed_url, session)
        if DOI is not None:
            row_record.append(DOI)
            return row_record
    return no_data_row

async def extract_from_db(url, session):
    async with semaphore:
        try:
            async with session.get(url, timeout=aiohttp.ClientTimeout(total=10)) as response:
                if response.status == 429:
                    retry_after = int(response.headers.get("Retry-After", 5))
                    print(f"Rate limited, sleeping {retry_after} seconds...")
                    await asyncio.sleep(retry_after)
                response.raise_for_status()  # raise error on bad status
                html_content = await response.text()
        except (aiohttp.ClientError, asyncio.TimeoutError) as e:
            print(f"Error fetching {url}: {e}")
            return None
    
        soup = BeautifulSoup(html_content, "html.parser")
        

        cell = soup.find("td", string="PubMed ID")
        pubmed_id = cell.find_next_sibling("td").get_text(strip=True)


        cell = soup.find("td", string = "Protein name")
        protein_name = cell.find_next_sibling("td").get_text(strip=True)

        cell = soup.find("td", string = "Comment")

        try:
            row = refoldDF[(refoldDF["pubmed_id"] == pubmed_id) & (refoldDF["name"] == protein_name)].iloc[0]
        
        except IndexError:
            return None


        refolding_method_columns =  ['dilution','dialysis','column_filtration','column_binding','high_pressure','other_method']

        # Apply only to subset of columns
        refolding_method = ', '.join(
        row[refolding_method_columns].index[row[refolding_method_columns].astype(bool)]
    )
        

        row_record = [row["name"], #Protein.Protein_Name
                      "", # Protein.Abbreviated_Name
                      "", # Protein.SCOP_Family
                      "", # Protein.Structure_Notes
                      "", # Protein.Organism
                      row["uniprot_id"], # Protein.UniProt_Accession
                      "", # Protein.SCOP_Unique_ID
                      "", # Protein.Structure_Solved
                      row["function"], # Protein.Class
                      "",# Protein.Molecularity
                      "", # Construct.Full_Length
                      row["domain"], # Construct.Domain
                      "", # Construct.Chimera
                      "", # Construct.Variants
                      len(row["aaseq"]) if pd.notna(row["aaseq"]) else "", # Construct.Chain_Length
                      ProteinAnalysis(row["aaseq"]).molecular_weight() if pd.notna(row["aaseq"]) else "", # Construct.Molecular_Weight
                      ProteinAnalysis(row["aaseq"]).isoelectric_point() if pd.notna(row["aaseq"]) else "", # Construct.Pi
                                            ProteinAnalysis(row["aaseq"]).molecular_weight() if pd.notna(row["aaseq"]) else "", # Construct.Molecular_Weight DUPLICATE ON ORIGNAL DB
                      "", # Construct.Disulphides
                      row["aaseq"], # Construct.Full_Sequence
                      row["comment"], # Construct.Notes
                      row["author"]+ row["date"] + " " + row["journal"], # Expression.Report
                      "", # Expression.Project_Aim
                      "", # Expression.Fusion
                      "", # Expression.Protein_Expression_and_Production
                      "", # Expression.Expression_Host
                      "", # Expression.Expression_Strain
                      "", # Expression.Expression_Temp
                      "", # Expression.Expression_Time
                      "", # Expression.Expression_Vector
                      "", # Expression.Expression_Protocol
                      "", # Expression.Method_of_Induction
                      "", # Expression.Cell_Density_at_Induction
                      "", # Expression.Cell_Disruption_Method
                      "", # Expression.Lytic_Agent
                      "", # Expression.Pre-Refolding_Purification
                      "", # Expression.Solubility
                      refolding_method, # Refolding.Refolding_Method
                      "", # Refolding.Wash_Buffer
                      "", # Refolding.Solubilization_Buffer
                      "", # Refolding.Refolding_Buffer
                      "", # Refolding.Pre-Refolding_Purification
                      "", # Refolding.Tag_Cleaved
                      row["ph"], # Refolding.Refolding_pH
                      row["temperature"], # Refolding.Refolding_Temperature
                      "", # Refolding.Protein_Concentration
                      "", # Refolding.Refolding_Time
                      "", # Refolding.Redox_Agent
                      "", # Refolding.Redox_Agent_Concentration
                      "", # Refolding.Refolding_Protocol
                      "", # Refolding.Refolding_Assay
                      "", # Refolding.Refolding_Chaperones
                      "", # Refolding.Refolding_Additives
                      "", # Refolding.Additives_Concentration
                      "", # Refolding.Refolding_Yield
                      "", # Refolding.Purity
                      "", # Refolding.Notes
        ]
        return row_record

async def main():

    async with aiohttp.ClientSession() as session:

        tasks = [process_url(url, session) for url in links]
        # Using tqdm with asyncio.gather
        links_values = []
        for f in tqdm_asyncio.as_completed(tasks, total=len(tasks), desc="Processing URLs"):
            result = await f
            links_values.append(result)
        return links_values

if __name__ == "__main__":
    start = datetime.now()
    loop = asyncio.get_event_loop()
    results = loop.run_until_complete(main())
    df = pd.DataFrame(results, columns=col_names)
    df = df.loc[:, ~df.columns.duplicated()] # Deal with duplicate Molecular_weight name in the data

    conn = sqlite3.connect('REFOLD.db')
    df.to_sql("FullDBTrial5",con=conn, index=False)