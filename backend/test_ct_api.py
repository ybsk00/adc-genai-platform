import asyncio
import aiohttp
import json

async def test_ct_api():
    # Search for a drug to get NCT IDs
    drug_name = "Trastuzumab deruxtecan"
    base_url = "https://clinicaltrials.gov/api/v2/studies"
    
    params = {
        "query.id": "NCT03248492",
        "pageSize": 1,
        "fields": "NCTId,BriefTitle,ResultsSection"
    }
    
    print(f"🔍 Searching ClinicalTrials.gov for {drug_name}...")
    
    async with aiohttp.ClientSession() as session:
        async with session.get(base_url, params=params) as response:
            if response.status != 200:
                print(f"❌ Error: {response.status}")
                text = await response.text()
                print(text)
                return

            data = await response.json()
            studies = data.get("studies", [])
            print(f"📚 Found {len(studies)} studies.")
            
            for study in studies:
                protocol = study.get("protocolSection", {})
                id_module = protocol.get("identificationModule", {})
                nct_id = id_module.get("nctId")
                title = id_module.get("briefTitle")
                
                print(f"\n📄 Study: {nct_id} - {title}")
                
                # Check Outcomes
                outcomes = protocol.get("outcomesModule", {}).get("primaryOutcomes", [])
                results = study.get("resultsSection", {})
                
                if results:
                    print("   ✅ Results Section Found!")
                    # Try to find outcome measures in results
                    outcome_measures = results.get("outcomeMeasuresModule", {}).get("outcomeMeasures", [])
                    for measure in outcome_measures[:3]: # Show first 3
                        print(f"      - Measure: {measure.get('title')}")
                        # Data is usually deep inside 'classes' -> 'categories' -> 'measurements'
                        # Just printing structure existence for now
                        print(f"        (Data structure present: {bool(measure.get('classes'))})")
                else:
                    print("   ⚠️ No Results Section (might be just protocol or no results posted yet)")
                    
                # Check Primary Outcomes (Descriptions)
                for outcome in outcomes[:2]:
                     print(f"      - Primary Outcome: {outcome.get('measure')}")

if __name__ == "__main__":
    asyncio.run(test_ct_api())
