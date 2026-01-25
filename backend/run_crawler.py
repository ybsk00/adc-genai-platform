import asyncio
import os
import sys
import logging
import argparse
from datetime import datetime
from dotenv import load_dotenv

# Add the backend directory to sys.path
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

# Load environment variables BEFORE importing any app modules
env_path = os.path.join(current_dir, ".env")
if os.path.exists(env_path):
    load_dotenv(env_path)
    print(f"✅ Loaded .env from {env_path}")
else:
    print(f"⚠️ .env file not found at {env_path}")

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("IsolatedCrawler")

async def main():
    parser = argparse.ArgumentParser(description="Run Crawler in an isolated process")
    parser.add_argument("--crawler", required=True, choices=["ambeed", "creative_biolabs"], help="Crawler type")
    parser.add_argument("--category", required=True, help="Category or Search term")
    parser.add_argument("--limit", type=int, default=10, help="Limit")
    parser.add_argument("--job_id", required=True, help="Job ID")
    parser.add_argument("--start_page", type=int, default=1, help="Start Page")
    
    args = parser.parse_args()
    
    logger.info(f"🚀 [Isolated Crawler] Starting {args.crawler} for {args.category} (Job: {args.job_id}, Page: {args.start_page})")
    
    try:
        if args.crawler == "ambeed":
            from app.services.ambeed_crawler import ambeed_crawler
            await ambeed_crawler.run(args.category, args.limit, args.job_id, args.start_page)
        elif args.crawler == "creative_biolabs":
            from app.services.creative_biolabs_crawler import creative_crawler
            # creative_biolabs doesn't have _run_internal yet, it uses run()
            # But run() in creative_biolabs DOES NOT spawn a subprocess, so it's safe to call here.
            await creative_crawler.run(args.category, args.limit, args.job_id)
            
    except Exception as e:
        logger.error(f"🔥 [Isolated Crawler] Failed: {e}", exc_info=True)
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
    # Force gRPC settings for this process to prevent fork issues
    os.environ["GRPC_ENABLE_FORK_SUPPORT"] = "false"
    os.environ["GRPC_POLL_STRATEGY"] = "epoll1"
    
    asyncio.run(main())