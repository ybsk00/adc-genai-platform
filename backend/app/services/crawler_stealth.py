import requests
from bs4 import BeautifulSoup
import time
import random
from fake_useragent import UserAgent
from supabase import create_client, Client
import os
from dotenv import load_dotenv
import logging

# Load environment variables
load_dotenv()

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class AmbeedStealthScraper:
    BASE_URL = "https://www.ambeed.com"
    
    def __init__(self):
        self.ua = UserAgent()
        self.session = requests.Session()
        
        # Initialize Supabase client
        url: str = os.environ.get("SUPABASE_URL")
        key: str = os.environ.get("SUPABASE_SERVICE_KEY")
        if not url or not key:
            raise ValueError("Supabase credentials not found in environment variables")
        self.supabase: Client = create_client(url, key)
        
    def get_headers(self, referer=None):
        headers = {
            'User-Agent': self.ua.random,
            'Accept-Language': 'en-US,en;q=0.9',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
            'Connection': 'keep-alive',
            'Upgrade-Insecure-Requests': '1',
        }
        if referer:
            headers['Referer'] = referer
        return headers

    def safe_request(self, url, referer=None):
        """
        Anti-blocking request wrapper with random delays and exponential backoff
        """
        retries = 3
        backoff = 30
        
        for i in range(retries):
            try:
                # 1. Human-like random delay (2-5 seconds)
                sleep_time = random.uniform(2, 5)
                time.sleep(sleep_time)
                
                logger.info(f"🕵️ Scraping: {url} (Delay: {sleep_time:.2f}s)")
                
                response = self.session.get(
                    url, 
                    headers=self.get_headers(referer), 
                    timeout=20
                )
                
                # 2. Handle blocking (403, 429, 503)
                if response.status_code in [403, 429, 503]:
                    logger.warning(f"⚠️ Blocked ({response.status_code}). Cooling down for {backoff}s...")
                    time.sleep(backoff)
                    backoff *= 2 # Exponential backoff
                    self.session = requests.Session() # Reset session (new identity)
                    continue
                
                response.raise_for_status()
                return response
            
            except Exception as e:
                logger.error(f"❌ Error fetching {url}: {str(e)}")
                time.sleep(5)
                
        return None

    def parse_detail_page(self, url, category, referer_url):
        """
        Parse product detail page and save to DB
        """
        res = self.safe_request(url, referer=referer_url)
        if not res:
            return
        
        soup = BeautifulSoup(res.text, 'html.parser')
        
        details = {
            'category': category,
            'product_url': url,
            'price_data': []
        }
        
        try:
            # Parse Product Name
            title_elem = soup.select_one('h1')
            if title_elem:
                details['product_name'] = title_elem.get_text(strip=True)

            # Parse Catalog Number (usually in breadcrumb or title)
            # This is site-specific, adjusting based on typical Ambeed structure
            # Assuming cat no is often prominent
            
            # Parse Specifications Table
            # Looking for table rows with specific labels
            for row in soup.select('table tr'):
                text = row.get_text(" ", strip=True)
                if "Catalog No" in text:
                    details['ambeed_cat_no'] = text.split(':')[-1].strip()
                elif "CAS No" in text:
                    details['cas_number'] = text.split(':')[-1].strip()
                elif "Molecular Formula" in text:
                    details['formula'] = text.split(':')[-1].strip()
                elif "Molecular Weight" in text:
                    details['molecular_weight'] = text.split(':')[-1].strip()
                elif "SMILES" in text: # Sometimes labeled as SMILES Code
                    details['smiles_code'] = text.split(':')[-1].strip()

            # Parse Pricing (if available in a table)
            # This requires inspecting the specific pricing table structure
            # Placeholder for pricing logic
            
            # Save to DB
            if 'ambeed_cat_no' in details:
                self.save_to_db(details)
            else:
                logger.warning(f"Skipping {url}: Could not find Catalog No")

        except Exception as e:
            logger.error(f"Error parsing {url}: {str(e)}")

    def save_to_db(self, data):
        """
        Upsert data into Supabase
        """
        try:
            # Check if exists
            existing = self.supabase.table('commercial_reagents')\
                .select('id')\
                .eq('ambeed_cat_no', data['ambeed_cat_no'])\
                .execute()
                
            if existing.data:
                # Update
                self.supabase.table('commercial_reagents')\
                    .update(data)\
                    .eq('ambeed_cat_no', data['ambeed_cat_no'])\
                    .execute()
                logger.info(f"Updated {data['ambeed_cat_no']}")
            else:
                # Insert
                self.supabase.table('commercial_reagents')\
                    .insert(data)\
                    .execute()
                logger.info(f"Inserted {data['ambeed_cat_no']}")
                
            # Trigger Embedding (Async in production, sync here for simplicity or via another service)
            # self.generate_embedding(data) 

        except Exception as e:
            logger.error(f"DB Error: {str(e)}")

    def run(self, category="Payload", max_pages=1):
        """
        Main execution loop
        """
        # Mapping categories to Ambeed search URLs (Placeholder URLs)
        category_urls = {
            "Payload": "https://www.ambeed.com/search?keyword=ADC+Payload",
            "Linker": "https://www.ambeed.com/search?keyword=ADC+Linker",
            "Conjugate": "https://www.ambeed.com/search?keyword=Antibody-Drug+Conjugate"
        }
        
        start_url = category_urls.get(category)
        if not start_url:
            logger.error(f"Invalid category: {category}")
            return

        logger.info(f"🚀 Starting crawler for {category}...")
        
        # 1. Visit List Page
        res = self.safe_request(start_url)
        if not res: return
        
        soup = BeautifulSoup(res.text, 'html.parser')
        
        # 2. Extract Product Links
        # Selector needs to be adjusted to actual Ambeed search results
        product_links = []
        for a in soup.select('.product-list a.title'): # Hypothetical selector
            href = a.get('href')
            if href:
                full_url = self.BASE_URL + href if href.startswith('/') else href
                product_links.append(full_url)
        
        logger.info(f"Found {len(product_links)} products on page 1")
        
        # 3. Visit Each Product
        for link in product_links:
            self.parse_detail_page(link, category, start_url)
            
        logger.info("✅ Crawling finished.")

if __name__ == "__main__":
    scraper = AmbeedStealthScraper()
    scraper.run()
