import asyncio
import os
import sys
from dotenv import load_dotenv

# Add backend directory to sys.path
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

load_dotenv(os.path.join(current_dir, ".env"))

from app.services.bulk_importer import BulkImporter

async def main():
    print("🚀 Starting ClinicalTrials.gov Sync...")
    importer = BulkImporter()
    # mode: 'daily' (updates) or 'full' (all)
    mode = "daily" 
    if len(sys.argv) > 1:
        mode = sys.argv[1]
    
    job_id = f"cmd_sync_clinical_{os.urandom(4).hex()}"
    print(f"Job ID: {job_id} | Mode: {mode}")
    
    await importer.run_import(job_id=job_id, max_studies=500, mode=mode)
    print("✅ ClinicalTrials.gov Sync Complete.")

if __name__ == "__main__":
    asyncio.run(main())
