import asyncio
import logging
import sys
import os

# Add backend directory to sys.path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from app.services.ambeed_crawler import ambeed_crawler

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def main():
    print("Starting Ambeed Crawler Debug...")
    
    # Test with ADC Toxins category, Page 1
    category = "ADC Toxins"
    url = ambeed_crawler.CATEGORIES[category]
    
    print(f"Crawling {category} from {url}...")
    
    # Run for just 1 page, limit 5 items
    count = await ambeed_crawler.crawl_category(
        category_name=category,
        base_url=url,
        limit=5,
        job_id=None, # No job ID to avoid DB updates for job status
        start_page=1,
        batch_size=2
    )
    
    print(f"Finished. Processed {count} items.")

if __name__ == "__main__":
    asyncio.run(main())
