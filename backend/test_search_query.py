import sys
import os
from datetime import datetime
from Bio import Entrez

# Add backend directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from app.core.config import settings

Entrez.email = settings.NCBI_EMAIL
if settings.NCBI_API_KEY:
    Entrez.api_key = settings.NCBI_API_KEY

SEARCH_QUERY = (
    '("Antibody-Drug Conjugate" OR "ADC" OR "Immunoconjugate" OR '
    '"Antibodies, Monoclonal, Conjugated"[MeSH] OR '
    '*-vedotin OR *-deruxtecan OR *-emtansine OR *-ozogamicin OR '
    '*-mafodotin OR *-tesirine) '
    'AND ("Clinical Trial" OR "Phase 1" OR "Phase 2" OR "Phase 3")'
)

print(f"Testing Query: {SEARCH_QUERY}")

try:
    handle = Entrez.esearch(
        db="pubmed",
        term=SEARCH_QUERY,
        retmax=5,
        sort="date",
        datetype="pdat",
        mindate=(datetime.now().year - 2),
        maxdate=datetime.now().year
    )
    record = Entrez.read(handle)
    handle.close()
    id_list = record.get("IdList", [])
    print(f"Found {len(id_list)} IDs: {id_list}")
except Exception as e:
    print(f"Error: {e}")
