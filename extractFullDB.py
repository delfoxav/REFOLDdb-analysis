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
import re

semaphore = asyncio.Semaphore(50)  # max 50 concurrent requests 
session = requests.Session() # Create a session so we don't create new ones later 


IDCONV_URL = " https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/"

MAILADRESS = "xavier@windyplace.ch" #requested for PMC ID Converter API https://pmc.ncbi.nlm.nih.gov/tools/id-converter-api/
TOOLNAME = "refoldbScrapper" #requested for PMC ID Converter API https://pmc.ncbi.nlm.nih.gov/tools/id-converter-api/
PMC_API_URL = f"https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/?tool={TOOLNAME}&email={MAILADRESS}"

BaseURL = "https://pford.info/refolddb/"
search_all = "./list.cgi?param=*&lang=EN" # An url to get all the entries from refoldDB website (the tsv is missing the pford links)
r = session.get("https://pford.info/refolddatabase/refoldingrecord/5/") # An example of an entry on pford (used to gather the keys of the new db)

soup = BeautifulSoup(r.text, "html.parser")
header = ""
col_names = []
data = []

for th in soup.find_all("th"):
    if th.get("class")[0] == "detail_header":
        header = th.find(string=True, recursive=False).strip()
    else:
        col_names.append(header + "." + th.text.replace(" ","_"))

col_names.append("DOI/PubMedID/PMCID")
print("Gathered all the col names")
check_paper_row = ["<UNKNOWN>"] * (len(col_names) - 1)
no_data_row = ["<UNKNOWN>"]*len(col_names)

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
    return ', '.join(row.index[row == 't'])

async def FindMonashUrl(url, session):
     """
     This asynchronous function is used to retrun the url to the Monash database from a refoldb entry.
     It keeps trying until  either an error 404 is reached or the url is found.

     If the refoldb entry has no url to Monash, it returns None
     """

     retry_after = 5
     while True: # retry forever
        async with semaphore:
            try:
                async with session.get(url, timeout=aiohttp.ClientTimeout(total=10)) as response:
                    
                    # if error 404 "not found" -> stop
                    if response.status == 404:
                        print(f"[404] Not found: {url}")
                        return None

                    # retry on any error
                    if response.status != 200:
                        retry_after = int(response.headers.get("Retry-After", 5))
                        print(f"Rate limited, sleeping {retry_after} seconds...")
                        await asyncio.sleep(retry_after)
                        continue # retry the request
                    
                    # if not error 429
                    response.raise_for_status()  # raise error on bad status

                    # in case of success
                    html_content = await response.text()

            except (aiohttp.ClientError, asyncio.TimeoutError) as e:
                print(f"Error fetching {url}: {e}, retrying in {retry_after} seconds... ")
                continue

            except Exception as e:
                # Unknown Error
                print(f"Fatal error for {url}: {e}")
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

async def FindPubMedID(url, session):
    """
     This asynchronous function is used to return the pubmed ID related to an article of the refoldb
     It keeps trying until  either an error 404 is reached or the ID is found.

     If the refoldb entry has no pubmed ID, it returns None
     """
    
    retry_after = 5
    while True: # retry forever
        async with semaphore:
            try:
                async with session.get(url, timeout=aiohttp.ClientTimeout(total=10)) as response:
                    
                    # if error 404 "not found" -> stop
                    if response.status == 404:
                        print(f"[404] Not found: {url}")
                        return None

                    # retry on any error
                    if response.status != 200:
                        retry_after = int(response.headers.get("Retry-After", 5))
                        print(f"Rate limited, sleeping {retry_after} seconds...")
                        await asyncio.sleep(retry_after)
                        continue # retry the request
                    
                    # if not error 429
                    response.raise_for_status()  # raise error on bad status

                    # in case of success
                    html_content = await response.text()

            except (aiohttp.ClientError, asyncio.TimeoutError) as e:
                print(f"Error fetching {url}: {e}, retrying in {retry_after} seconds... ")
                continue

            except Exception as e:
                # Unknown Error
                print(f"Fatal error for {url}: {e}")
                return None

        soup = BeautifulSoup(html_content, "html.parser")

        # Find all the table rows
        rows = soup.find_all("tr")

        for row in rows:
            # checks if one of the column is "PubMed ID"
            if row.find("td", string="PubMed ID") != None:
                if row.find("a") != None:
                    return row.find("a").text
                else:
                    return None
        
async def FindDOI_Batch(pmids, session, batch_size=200):
    """
    Batch DOI lookup using NCBI ID Converter API.
    Input:  list of PMIDs
    Output: dict { pmid: doi or None }
    """

    results = {str(pmid): None for pmid in pmids}  # fill defaults

    # Process in chunks (NCBI handles long lists, but chunking is safer)
    chunks = [
        pmids[i:i + batch_size]
        for i in range(0, len(pmids), batch_size)
    ]

    retry_after = 5

    for chunk in chunks:
        params = {
            "ids": ",".join(map(str, chunk)),
            "format": "json"
        }

        while True:
            try:
                async with session.get(IDCONV_URL, params=params,
                                       timeout=aiohttp.ClientTimeout(total=10)) as response:

                    if response.status != 200:
                        retry_after = int(response.headers.get("Retry-After", retry_after))
                        print(f"Batch rate-limited ({response.status}), sleeping {retry_after}s...")
                        await asyncio.sleep(retry_after)
                        continue

                    data = await response.json()

            except (aiohttp.ClientError, asyncio.TimeoutError) as e:
                print(f"Network error during batch DOI lookup: {e}; retrying in {retry_after}s...")
                await asyncio.sleep(retry_after)
                continue

            except Exception as e:
                print(f"Fatal batch DOI error: {e}")
                return results

            break  # success

        # ---------- Parse the batch ----------
        for record in data.get("records", []):
            pmid = record.get("pmid")
            doi = record.get("doi")
            pmcid = record.get("pmcid")

            if pmid is not None and pmcid is not None:
                results[str(pmid)] = pmcid
            elif pmid is not None and doi is not None:
                results[str(pmid)] = doi

    return results

async def gatherrecord(url, session):
    """
    This asynchronous function is used to return all the values of a single record from the monash DB
    It keeps trying until reaching an error 404.
    """
    retry_after = 5 # default value
    while True: # retry forever
        async with semaphore:
            try:
                async with session.get(url, timeout=aiohttp.ClientTimeout(total=10)) as response:
                    
                    # if error 404 "not found" -> stop
                    if response.status == 404:
                        print(f"[404] Not found: {url}")
                        return None

                    # retry on any error
                    if response.status != 200:
                        retry_after = int(response.headers.get("Retry-After", 5))
                        print(f"Rate limited, sleeping {retry_after} seconds...")
                        await asyncio.sleep(retry_after)
                        continue # retry the request
                    
                    # if not error 429
                    response.raise_for_status()  # raise error on bad status

                    # in case of success
                    html_content = await response.text()

            except (aiohttp.ClientError, asyncio.TimeoutError) as e:
                print(f"Error fetching {url}: {e}, retrying in {retry_after} seconds... ")
                continue

            except Exception as e:
                # Unknown Error
                print(f"Fatal error for {url}: {e}")
                return None
        soup = BeautifulSoup(html_content, "html.parser")

        data = []
        for row in soup.find_all('tr'):
                for cell in row.find_all("td"):
                    data.append(cell.text)
        
        return data

async def extract_from_db(url, session):
    retry_after = 5
    while True: # retry forever
        async with semaphore:
            try:
                async with session.get(url, timeout=aiohttp.ClientTimeout(total=10)) as response:
                    
                    # if error 404 "not found" -> stop
                    if response.status == 404:
                        print(f"[404] Not found: {url}")
                        return None

                    # retry on any error
                    if response.status != 200:
                        retry_after = int(response.headers.get("Retry-After", 5))
                        print(f"Rate limited, sleeping {retry_after} seconds...")
                        await asyncio.sleep(retry_after)
                        continue # retry the request
                    
                    # if not error 429
                    response.raise_for_status()  # raise error on bad status

                    # in case of success
                    html_content = await response.text()

            except (aiohttp.ClientError, asyncio.TimeoutError) as e:
                print(f"Error fetching {url}: {e}, retrying in {retry_after} seconds... ")
                continue

            except Exception as e:
                # Unknown Error
                print(f"Fatal error for {url}: {e}")
                return check_paper_row.copy()
    
        soup = BeautifulSoup(html_content, "html.parser")
        
        try:
            cell = soup.find("td", string="PubMed ID")
            pubmed_id = cell.find_next_sibling("td").get_text(strip=True)
        except AttributeError:
            print(f"No pubMed Id for {url}")
            return check_paper_row.copy()


        try:
            cell = soup.find("td", string = "Protein name")
            protein_name = cell.find_next_sibling("td").get_text(strip=True)
        except AttributeError:
            print(f"No protein Name for {url}")
            return check_paper_row.copy()

        try:
            row = refoldDF[(refoldDF["pubmed_id"] == pubmed_id) & (refoldDF["name"] == protein_name)].iloc[0]
        
        except IndexError:
            return check_paper_row.copy()

        refolding_method_columns = ['dilution', 'dialysis', 'column_filtration', 'column_binding', 'high_pressure', 'other_method']
        refolding_method = true_columns(row[refolding_method_columns])
        chain_length = len(row["aaseq"]) if pd.notna(row["aaseq"]) else "<UNKNOWN>"

    # Calculate Construct.Molecular_Weight and Construct.Pi with error handling
        if pd.notna(row["aaseq"]):
            try:
                pa = ProteinAnalysis(row["aaseq"])
                molecular_weight = pa.molecular_weight()
            except Exception:
                molecular_weight = "<UNKNOWN>"
            try:
                pi = pa.isoelectric_point()
            except Exception:
                pi = "<UNKNOWN>"
            # Duplicate for original DB
            try:
                molecular_weight_dup = pa.molecular_weight()
            except Exception:
                molecular_weight_dup = "<UNKNOWN>"
        else:
            molecular_weight = "<UNKNOWN>"
            pi = "<UNKNOWN>"
            molecular_weight_dup = "<UNKNOWN>"

        row_record = [
            row["name"],  # Protein.Protein_Name
            "<UNKNOWN>",  # Protein.Abbreviated_Name
            "<UNKNOWN>",  # Protein.SCOP_Family
            "<UNKNOWN>",  # Protein.Structure_Notes
            "<UNKNOWN>",  # Protein.Organism
            row["uniprot_id"],  # Protein.UniProt_Accession
            "<UNKNOWN>",  # Protein.SCOP_Unique_ID
            "<UNKNOWN>",  # Protein.Structure_Solved
            row["function"],  # Protein.Class
            "<UNKNOWN>",  # Protein.Molecularity
            "<UNKNOWN>",  # Construct.Full_Length
            row["domain"],  # Construct.Domain
            "<UNKNOWN>",  # Construct.Chimera
            "<UNKNOWN>",  # Construct.Variants
            chain_length,  # Construct.Chain_Length
            molecular_weight,  # Construct.Molecular_Weight
            pi,  # Construct.Pi
            molecular_weight_dup,  # Construct.Molecular_Weight DUPLICATE ON ORIGINAL DB
            "<UNKNOWN>",  # Construct.Disulphides
            row["aaseq"],  # Construct.Full_Sequence
            row["comment"],  # Construct.Notes
            row["author"] + row["date"] + " " + row["journal"],  # Expression.Report
            "<UNKNOWN>",  # Expression.Project_Aim
            "<UNKNOWN>",  # Expression.Fusion
            "<UNKNOWN>",  # Expression.Protein_Expression_and_Production
            "<UNKNOWN>",  # Expression.Expression_Host
            "<UNKNOWN>",  # Expression.Expression_Strain
            "<UNKNOWN>",  # Expression.Expression_Temp
            "<UNKNOWN>",  # Expression.Expression_Time
            "<UNKNOWN>",  # Expression.Expression_Vector
            "<UNKNOWN>",  # Expression.Expression_Protocol
            "<UNKNOWN>",  # Expression.Method_of_Induction
            "<UNKNOWN>",  # Expression.Cell_Density_at_Induction
            "<UNKNOWN>",  # Expression.Cell_Disruption_Method
            "<UNKNOWN>",  # Expression.Lytic_Agent
            "<UNKNOWN>",  # Expression.Pre-Refolding_Purification
            "<UNKNOWN>",  # Expression.Solubility
            refolding_method,  # Refolding.Refolding_Method
            "<UNKNOWN>",  # Refolding.Wash_Buffer
            "<UNKNOWN>",  # Refolding.Solubilization_Buffer
            "<UNKNOWN>",  # Refolding.Refolding_Buffer
            "<UNKNOWN>",  # Refolding.Pre-Refolding_Purification
            "<UNKNOWN>",  # Refolding.Tag_Cleaved
            row["ph"],  # Refolding.Refolding_pH
            row["temperature"],  # Refolding.Refolding_Temperature
            "<UNKNOWN>",  # Refolding.Protein_Concentration
            "<UNKNOWN>",  # Refolding.Refolding_Time
            "<UNKNOWN>",  # Refolding.Redox_Agent
            "<UNKNOWN>",  # Refolding.Redox_Agent_Concentration
            "<UNKNOWN>",  # Refolding.Refolding_Protocol
            "<UNKNOWN>",  # Refolding.Refolding_Assay
            "<UNKNOWN>",  # Refolding.Refolding_Chaperones
            "<UNKNOWN>",  # Refolding.Refolding_Additives
            "<UNKNOWN>",  # Refolding.Additives_Concentration
            "<UNKNOWN>",  # Refolding.Refolding_Yield
            "<UNKNOWN>",  # Refolding.Purity
            "<UNKNOWN>",  # Refolding.Notes
        ]
        return row_record



async def process_url(url, session, doi_map):
    monash_url = await FindMonashUrl(url, session) # first we check if the entry is in the monash db as well
    pubmed_id = await FindPubMedID(url, session) # We check if we have the PubMedID
    DOI = doi_map.get(str(pubmed_id))
    

    if monash_url is not None:
        row_record = await gatherrecord(monash_url, session) # if we have a monash url we can simply use the record
        if DOI is not None: # If we have a DOI we can add it to our record (giving us a ref for the info of that record)
            row_record.append(DOI)
        elif pubmed_id is not None: # otherwise we might have a pubmed_url that we can give to our record as a reference for info
            row_record.append(pubmed_id)
        else: # Lastly we might not have any. If so we can't really get a reference for the info like that
            row_record.append("<UNKNOWN>")
        return row_record
    
    else: # if  we don't have a monash url we will extract the info from refoldb as much as we can
        row_record = await extract_from_db(url, session)
        if DOI is not None: # If we have a DOI we can add it to our record (giving us a ref for the info of that record)
            row_record.append(DOI)
        elif pubmed_id is not None: # otherwise we might have a pubmed_url that we can give to our record as a reference for info
            row_record.append(pubmed_id)
        else: # Lastly we might not have any. If so we can't really get a reference for the info like that
            row_record.append("<UNKNOWN>")
        return row_record


async def main():

    async with aiohttp.ClientSession() as session:
        
        # First we gather all the DOIs with pmc api
        pmid_tasks = [FindPubMedID(url, session) for url in links]
        pmids_raw = []
        for f in tqdm_asyncio.as_completed(pmid_tasks, total=len(pmid_tasks), desc="Getting PubMedIds"):
            pmids_raw.append(await f)
        
        # Deduplicate + drop Nones
        pmids = sorted({pmid for pmid in pmids_raw if pmid is not None})

        doi_map = await FindDOI_Batch(pmids, session)
        print("built doi_map")

        tasks = [process_url(url, session, doi_map) for url in links]
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
    df.to_sql("FullDB",con=conn, index=False)