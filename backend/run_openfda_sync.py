import asyncio
import os
import sys
from dotenv import load_dotenv

# Add backend directory to sys.path
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

load_dotenv(os.path.join(current_dir, ".env"))

from app.services.openfda_service import openfda_service

async def main():
    print("🚀 Starting OpenFDA Sync...")
    # mode: 'daily' (new approvals) or 'full' (all)
    mode = "daily"
    limit = 100
    
    if len(sys.argv) > 1:
        mode = sys.argv[1]
    if len(sys.argv) > 2:
        limit = int(sys.argv[2])
    
    job_id = f"cmd_sync_openfda_{os.urandom(4).hex()}"
    print(f"Job ID: {job_id} | Mode: {mode} | Limit: {limit}")
    
    await openfda_service.sync_to_db(job_id=job_id, mode=mode, limit=limit)
    print("✅ OpenFDA Sync Complete.")

if __name__ == "__main__":
    asyncio.run(main())
