import asyncio
import sys
import os

# Add backend directory to sys.path to allow imports from app
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from app.api.scheduler import process_clinical_trials_data, sync_jobs, fetch_clinical_trials
from app.core.supabase import supabase

async def test_sync_logic():
    print("--- Testing Scheduler Logic ---")
    job_id = "test_job_1"
    
    # Initialize job in sync_jobs as the API would
    sync_jobs[job_id] = {
        "status": "queued",
        "source": "clinical_trials",
        "records_found": 0,
        "records_drafted": 0,
        "started_at": "now",
        "errors": []
    }

    print(f"Starting job {job_id}...")
    try:
        await process_clinical_trials_data(job_id)
        
        job_result = sync_jobs[job_id]
        print(f"Job Finished. Status: {job_result['status']}")
        print(f"Records Found: {job_result['records_found']}")
        print(f"Records Drafted: {job_result['records_drafted']}")
        
        if job_result['errors']:
            print("Errors encountered:")
            for err in job_result['errors']:
                print(f" - {err}")
        
        # Verify Supabase
        print("\nVerifying Supabase Data...")
        response = supabase.table("golden_set_library").select("count", count="exact").execute()
        print(f"Total records in golden_set_library: {response.count}")

    except Exception as e:
        print(f"CRITICAL ERROR: {e}")
        if hasattr(e, 'response'):
            print(f"Response Status: {e.response.status_code}")
            print(f"Response Text: {e.response.text}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    asyncio.run(test_sync_logic())
