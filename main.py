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
    monash_url = await FindMonashUrl(url, session)
    if monash_url is not None:
        row_record = await gatherrecord(monash_url, session)
        row_record.append("Not needed")
        return row_record
    pubmed_url = await FindPubMedURL(url, session)
    if pubmed_url is not None:
        row = check_paper_row.copy()
        DOI = await FindDOI(pubmed_url, session)
        if DOI is not None:
            row.append(DOI)
            return row
    return no_data_row

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
    df.to_sql("FullDBTrial",con=conn, index=False)