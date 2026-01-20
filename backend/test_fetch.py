import asyncio
import httpx
from Bio import Entrez
import os

# Entrez 설정
Entrez.email = "test@example.com"

async def test_clinical_trials():
    print("\n--- Testing ClinicalTrials.gov API ---")
    query = "ADC OR Antibody-Drug Conjugate"
    url = "https://clinicaltrials.gov/api/v2/studies"
    params = {
        "query.term": query,
        "pageSize": 5,
        "format": "json"
    }
    try:
        async with httpx.AsyncClient() as client:
            print(f"Fetching from {url} with params: {params}")
            response = await client.get(url, params=params, timeout=30.0)
            print(f"Status Code: {response.status_code}")
            if response.status_code == 200:
                data = response.json()
                studies = data.get("studies", [])
                print(f"Found {len(studies)} studies.")
                if studies:
                    print(f"First study title: {studies[0].get('protocolSection', {}).get('identificationModule', {}).get('officialTitle')}")
            else:
                print(f"Error: {response.text}")
    except Exception as e:
        print(f"Exception: {e}")

def test_pubmed():
    print("\n--- Testing PubMed API ---")
    term = "Antibody-Drug Conjugate"
    try:
        print(f"Searching PubMed for: {term}")
        handle = Entrez.esearch(db="pubmed", term=term, retmax=5)
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]
        print(f"Found {len(id_list)} IDs: {id_list}")
        
        if id_list:
            print("Fetching details for first ID...")
            handle = Entrez.efetch(db="pubmed", id=id_list[0], retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            article = records['PubmedArticle'][0]
            title = article['MedlineCitation']['Article']['ArticleTitle']
            print(f"First article title: {title}")
    except Exception as e:
        print(f"Exception: {e}")

if __name__ == "__main__":
    asyncio.run(test_clinical_trials())
    test_pubmed()
