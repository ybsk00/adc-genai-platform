import asyncio
import os
import sys
import logging
import argparse
from datetime import datetime

# Add the backend directory to sys.path so we can import app modules
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

async def main():
    parser = argparse.ArgumentParser(description="Run Ambeed Crawler in an isolated process")
    parser.add_argument("--job_id", required=True, help="Job ID")
    parser.add_argument("--search_term", required=True, help="Search term")
    parser.add_argument("--limit", type=int, default=20, help="Limit")
    
    args = parser.parse_args()
    
    logger.info(f"🚀 [Isolated Crawler] Starting for Job: {args.job_id}")
    
    try:
        # Import here to avoid early initialization issues
        from app.services.ambeed_crawler import ambeed_crawler
        
        # We need to bypass the subprocess spawning logic in run() to avoid infinite recursion
        # We'll create a new method _run_internal in ambeed_crawler.py or just call the logic if we can.
        # But wait, I haven't modified ambeed_crawler.py yet.
        # I will modify ambeed_crawler.py to have _run_internal and run.
        # run() will spawn subprocess. _run_internal() will do the work.
        # This script will call _run_internal().
        
        await ambeed_crawler._run_internal(args.search_term, args.limit, args.job_id)
        
    except Exception as e:
        logger.error(f"🔥 [Isolated Crawler] Failed: {e}", exc_info=True)
        # We should probably update the job status to failed here if ambeed_crawler didn't
        try:
            from app.core.supabase import supabase
            supabase.table("sync_jobs").update({
                "status": "failed",
                "errors": [str(e)]
            }).eq("id", args.job_id).execute()
        except:
            pass
        sys.exit(1)

if __name__ == "__main__":
    # Force gRPC settings for this process
    os.environ["GRPC_ENABLE_FORK_SUPPORT"] = "false"
    os.environ["GRPC_POLL_STRATEGY"] = "epoll1"
    
    asyncio.run(main())
