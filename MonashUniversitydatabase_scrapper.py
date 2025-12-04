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

urls = []
session = requests.Session()

# first we get the col_names from a single record
r = requests.get("https://pford.info/refolddatabase/refoldingrecord/5/")

soup = BeautifulSoup(r.text, "html.parser")
header = ""
col_names = []
data = []

for th in soup.find_all("th"):
    if th.get("class")[0] == "detail_header":
        header = th.find(string=True, recursive=False).strip()
    else:
        col_names.append(header + "." + th.text.replace(" ","_"))

print("Gathered all the col names")

# We get a list of pages urls

r = session.get("https://pford.info/refolddatabase/refoldingrecord")
soup = BeautifulSoup(r.text, "html.parser")
while True:
    next_link = soup.find("a", string="Next")
    if not next_link:
        break
    url = "https://pford.info/refolddatabase/refoldingrecord" + next_link.get("href")
    urls.append(url)
    
    r = session.get(url)
    soup = BeautifulSoup(r.text, "html.parser")
    

print("Gathered all the urls")

async def gatherrecord(url):
    """
    Function to gather all the value of a single record
    """

    async with aiohttp.ClientSession() as session:
        async with session.get(url) as response:
            html_content = await response.text()

    soup = BeautifulSoup(html_content, "html.parser")

    data = []
    for row in soup.find_all('tr'):
            for cell in row.find_all("td"):
                data.append(cell.text)
    
    return data


async def analyze_page(url):
    """
    Function to get all the records of a page
    """

    async with aiohttp.ClientSession() as session:
        async with session.get(url) as response:
            html_content = await response.text()

    soup = BeautifulSoup(html_content, "html.parser")

    page_data = []
    for row in tqdm(soup.find_all('tr'), unit="record"):
        link = row.find("a")
        if link: # Deal with the headers 
            url = "https://pford.info/refolddatabase/refoldingrecord/" + link.get("href")
            row_data =  await gatherrecord(url)
            page_data.append(row_data)
                

    return page_data

async def main():
    tasks = [analyze_page(url) for url in urls]
    pages = await asyncio.gather(*tasks)   # pages: list of pages
    # flatten pages -> list of records
    records = [record for page in pages for record in page]
    return records



if __name__ == "__main__":
    start = datetime.now()
    loop = asyncio.get_event_loop()
    results = loop.run_until_complete(main())
    print("Total Time: ", datetime.now()-start)
    
    df = pd.DataFrame(results, columns=col_names)
    df = df.loc[:, ~df.columns.duplicated()] # Deal with duplicate Molecular_weight name in the data
    conn = sqlite3.connect('REFOLD.db')
    df.to_sql("main",con=conn, index=False)