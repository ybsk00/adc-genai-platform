import asyncio
import argparse
import json
import logging
import os
import sys
import re
import aiohttp
from typing import List, Dict, Any, Optional
from datetime import datetime

# Add backend directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from app.core.supabase import supabase
from app.core.config import settings
from app.services.pubmed_knowledge_service import pubmed_knowledge_service
import google.generativeai as genai
from json_repair import repair_json

# Configure Logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

class GoldenSetEnricher:
    """
    ADC Golden Set Enricher & Expander
    Strategy:
    1. Triple-Layered Search Query (MeSH + Suffix + Keywords)
    2. 4-Point AI Extraction (Target, Antibody, Linker, Payload)
    3. ClinicalTrials.gov v2 Integration (Metric Booster)
    4. Strict Upsert Logic (Update Nulls Only vs Insert New)
    5. Easy Win SMILES
    """

    # 1. Triple-Layered Search Query
    SEARCH_QUERY = (
        '("Antibody-Drug Conjugate" OR "ADC" OR "Immunoconjugate" OR '
        '"Antibodies, Monoclonal, Conjugated"[MeSH] OR '
        '*-vedotin OR *-deruxtecan OR *-emtansine OR *-ozogamicin OR '
        '*-mafodotin OR *-tesirine) '
        'AND ("Clinical Trial" OR "Phase 1" OR "Phase 2" OR "Phase 3")'
    )

    def __init__(self):
        genai.configure(api_key=settings.GOOGLE_API_KEY)
        self.model = genai.GenerativeModel('gemini-2.0-flash')
        self.processed_count = 0
        self.updated_count = 0
        self.new_count = 0

    async def run_discovery(self, max_results: int = 20, dry_run: bool = False):
        """Run discovery mode using the Triple-Layered Query"""
        logger.info(f"🚀 Starting Discovery Mode with Query: {self.SEARCH_QUERY}")
        
        from Bio import Entrez
        Entrez.email = settings.NCBI_EMAIL
        if settings.NCBI_API_KEY:
            Entrez.api_key = settings.NCBI_API_KEY

        try:
            handle = Entrez.esearch(
                db="pubmed",
                term=self.SEARCH_QUERY,
                retmax=max_results,
                sort="date", # Get latest first
                datetype="pdat",
                mindate=(datetime.now().year - 2), # Last 2 years for "New" stuff
                maxdate=datetime.now().year
            )
            record = Entrez.read(handle)
            handle.close()
            id_list = record.get("IdList", [])
            
            logger.info(f"📚 Found {len(id_list)} articles.")
            
            if not id_list:
                return

            # Fetch details
            articles = await pubmed_knowledge_service._fetch_details(id_list)
            
            for article in articles:
                await self.process_article(article, dry_run)
                
        except Exception as e:
            logger.error(f"Search failed: {e}")

    async def process_article(self, article: Dict, dry_run: bool):
        """Process a single article: AI Extract -> CT.gov -> Upsert"""
        title = article['title']
        abstract = article['abstract']
        pmid = article['pmid']
        
        logger.info(f"📄 Processing: {title[:50]}... (PMID: {pmid})")
        
        # 2. AI Extraction (4-Point Focus)
        extracted_data = await self.extract_adc_info(title, abstract)
        
        if not extracted_data:
            # logger.warning("   ⚠️ AI Extraction failed or irrelevant.") # Already logged inside
            return

        # 3. ClinicalTrials.gov v2 Integration (Metric Booster)
        # Try to find NCT ID in abstract
        nct_match = re.search(r'NCT\d{8}', abstract)
        if nct_match:
            nct_id = nct_match.group(0)
            logger.info(f"      🔗 Found NCT ID: {nct_id}. Fetching official metrics...")
            ct_metrics = await self.fetch_clinical_metrics(nct_id)
            if ct_metrics:
                logger.info(f"      📊 Overwriting AI metrics with Official Data: {ct_metrics}")
                extracted_data.update(ct_metrics)

        if dry_run:
            logger.info(f"   [DRY RUN] Extracted: {json.dumps(extracted_data, indent=2)}")
            return

        # 4. Upsert Logic
        await self.upsert_record(extracted_data)

    async def fetch_clinical_metrics(self, nct_id: str) -> Optional[Dict]:
        """
        Fetch official metrics from ClinicalTrials.gov API v2
        """
        url = f"https://clinicaltrials.gov/api/v2/studies/{nct_id}"
        params = {"fields": "ResultsSection"}
        
        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params) as response:
                    if response.status != 200: return None
                    
                    data = await response.json()
                    results = data.get("resultsSection", {})
                    if not results: return None
                    
                    measures = results.get("outcomeMeasuresModule", {}).get("outcomeMeasures", [])
                    metrics = {}
                    
                    for m in measures:
                        title = m.get("title", "").lower()
                        try:
                            # Navigate deep structure: classes -> categories -> measurements -> value
                            # This path is common for simple numeric outcomes
                            value_str = m.get("classes", [])[0].get("categories", [])[0].get("measurements", [])[0].get("value")
                            if not value_str: continue
                            
                            # Clean value (remove confidence intervals etc)
                            val = float(re.split(r'\s|\[|\(', str(value_str))[0])
                            
                            if "response" in title and "rate" in title:
                                metrics["orr_pct"] = val
                            elif "progression" in title and "free" in title:
                                metrics["pfs_months"] = val
                            elif "overall" in title and "survival" in title:
                                metrics["os_months"] = val
                        except:
                            continue
                            
                    if metrics:
                        metrics["outcome_type"] = "Success" # If results are posted, it's likely a completed/successful trial phase
                        return metrics
                    return None
        except Exception as e:
            logger.error(f"      ❌ CT.gov API Error: {e}")
            return None

    async def extract_adc_info(self, title: str, abstract: str) -> Optional[Dict]:
        """
        Gemini 2.0 Flash Prompt for 4-Point Extraction & Metrics
        """
        prompt = f"""
        Analyze this ADC clinical trial abstract.
        
        Title: {title}
        Abstract: {abstract}
        
        Task: Extract structured ADC data. Focus on the "4-Point Structure" and "Clinical Metrics".
        
        Rules:
        1. **Identity**: Extract the Drug Name (e.g., "Trastuzumab deruxtecan", "DS-8201").
        2. **4-Point Structure**: Split into Target, Antibody, Linker, Payload.
           - Target: e.g., "HER2"
           - Antibody: e.g., "Trastuzumab" (or "IgG1")
           - Linker: e.g., "Tetrapeptide" or "Cleavable"
           - Payload: e.g., "DXd" or "MMAE"
           - *IMPORTANT*: If a component is implied (e.g., "vedotin" implies MMAE payload and cleavable linker), EXTRACT IT explicitly. Do not leave null if it can be inferred.
        3. **Metrics**: Extract NUMERIC values only.
           - ORR (Objective Response Rate) in % (e.g., extract 43.5 from "43.5%" or "ORR was 43.5%"). If a range is given, take the median or best value.
           - PFS (Progression Free Survival) in months.
           - OS (Overall Survival) in months.
           - DAR (Drug Antibody Ratio).
           - *CRITICAL*: Look closely for these numbers. They are the most important part.
        4. **Outcome**: Determine if the trial was a "Success" or "Failure".
           - Success: Met endpoints, approved, promising results, "manageable safety profile", "significant improvement".
           - Failure: Terminated, did not meet endpoints, high toxicity, "limited efficacy".
           - *DEFAULT*: If unsure, use "Unknown" but try to infer from the tone (positive=Success, negative=Failure).
        5. **SMILES**: ONLY if explicitly mentioned or commonly known (Easy Win). Otherwise null.
        
        Output JSON:
        {{
            "drug_name": "string",
            "target": "string",
            "antibody": "string",
            "linker": "string",
            "payload": "string",
            "orr_pct": float or null,
            "pfs_months": float or null,
            "os_months": float or null,
            "dar": float or null,
            "outcome_type": "Success" | "Failure" | "Unknown",
            "smiles": "string" or null
        }}
        """
        
        try:
            response = await self.model.generate_content_async(prompt)
            text = response.text
            
            # Repair and Parse JSON
            json_str = repair_json(text)
            data = json.loads(json_str)
            
            # 1. Handle List vs Dict (Bug Fix)
            if isinstance(data, list):
                if not data: return None
                data = data[0] # Take the first item if it's a list
            
            if not isinstance(data, dict):
                logger.warning("      ⚠️ AI returned non-dict data.")
                return None

            # 2. Filtering Logic (Irrelevant / Protocol)
            # Check for "Study Protocol" or irrelevant content
            if "protocol" in title.lower() or "design" in title.lower() and "rationale" in title.lower():
                logger.info("      ⏩ Skipping: Likely a Study Protocol (No results yet).")
                return None
            
            # Basic Validation: Must have a drug name and at least one ADC component
            if not data.get("drug_name") or not (data.get("target") or data.get("payload")):
                logger.warning("      ⚠️ Missing core ADC data (Drug Name/Target/Payload).")
                return None
                
            return data
        except Exception as e:
            logger.error(f"   ❌ AI Error: {e}")
            return None

    async def upsert_record(self, data: Dict):
        """
        Strict Upsert Logic:
        - Match by Name (fuzzy or exact)
        - If Exists: Update ONLY Null fields
        - If New: Insert Full Record
        """
        drug_name = data['drug_name']
        
        # Check existence
        # Using a simple ILIKE for now. In production, might need more sophisticated matching.
        res = supabase.table("golden_set_library").select("*").ilike("name", drug_name).execute()
        
        if res.data:
            # === UPDATE EXISTING ===
            existing = res.data[0]
            record_id = existing['id']
            logger.info(f"   🔄 Found existing record: {existing['name']} (ID: {record_id})")
            
            updates = {}
            properties = existing.get("properties") or {}
            
            # Map extracted fields to DB columns/properties
            # Columns: outcome_type, dar, orr_pct, os_months, pfs_months, target_1, antibody_format, linker_type, payload_smiles
            
            # Helper to update if null
            def update_if_null(db_key, new_val, is_prop=False):
                if new_val is None: return
                
                current_val = existing.get(db_key) if not is_prop else properties.get(db_key)
                
                if current_val is None or current_val == "":
                    if is_prop:
                        properties[db_key] = new_val
                        # Mark properties as modified
                        updates["properties"] = properties 
                    else:
                        updates[db_key] = new_val
                    logger.info(f"      + Updating {db_key}: {new_val}")

            # Apply updates
            update_if_null("target_1", data.get("target"))
            update_if_null("antibody_format", data.get("antibody"))
            update_if_null("linker_type", data.get("linker"))
            update_if_null("payload_smiles", data.get("smiles")) # Easy Win SMILES
            
            update_if_null("orr_pct", data.get("orr_pct"))
            update_if_null("pfs_months", data.get("pfs_months"))
            update_if_null("os_months", data.get("os_months"))
            update_if_null("dar", data.get("dar"))
            update_if_null("outcome_type", data.get("outcome_type"))
            
            # Payload Name (Property)
            update_if_null("payload", data.get("payload"), is_prop=True)

            if updates:
                supabase.table("golden_set_library").update(updates).eq("id", record_id).execute()
                self.updated_count += 1
                logger.info("      ✅ Updates applied.")
            else:
                logger.info("      Example: No null fields to update.")

        else:
            # === INSERT NEW ===
            logger.info(f"   ✨ New Drug Detected: {drug_name}")
            
            new_record = {
                "name": drug_name,
                "status": "draft", # Mark as draft for review
                "enrichment_source": "GoldenSetEnricher",
                "target_1": data.get("target"),
                "antibody_format": data.get("antibody"),
                "linker_type": data.get("linker"),
                "payload_smiles": data.get("smiles"),
                "orr_pct": data.get("orr_pct"),
                "pfs_months": data.get("pfs_months"),
                "os_months": data.get("os_months"),
                "dar": data.get("dar"),
                "outcome_type": data.get("outcome_type"),
                "properties": {
                    "payload": data.get("payload"),
                    "ai_extracted": True
                }
            }
            
            supabase.table("golden_set_library").insert(new_record).execute()
            self.new_count += 1
            logger.info("      ✅ New record inserted.")

async def main():
    parser = argparse.ArgumentParser(description="Golden Set Enricher")
    parser.add_argument("--mode", choices=["discovery", "target"], default="discovery", help="Mode: discovery (all) or target (specific drug)")
    parser.add_argument("--target", help="Target drug name for target mode")
    parser.add_argument("--dry-run", action="store_true", help="Dry run (no DB changes)")
    parser.add_argument("--limit", type=int, default=10, help="Max articles to process")
    
    args = parser.parse_args()
    
    enricher = GoldenSetEnricher()
    
    if args.mode == "discovery":
        await enricher.run_discovery(max_results=args.limit, dry_run=args.dry_run)
    elif args.mode == "target":
        # TODO: Implement target specific search if needed, 
        # for now discovery covers the general enrichment case.
        logger.info("Target mode not fully implemented yet, use discovery.")
        
    logger.info(f"\n🏁 Finished. Updated: {enricher.updated_count}, New: {enricher.new_count}")

if __name__ == "__main__":
    asyncio.run(main())
