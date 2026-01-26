import asyncio
import logging
import json
import os
import re
from datetime import datetime
from dotenv import load_dotenv
from playwright.async_api import async_playwright

load_dotenv()
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("TestDirectSave")

# --- Supabase Init ---
from supabase import create_client
url = os.getenv("SUPABASE_URL")
key = os.getenv("SUPABASE_SERVICE_KEY")
supabase = create_client(url, key)

async def test_single_product():
    # 사장님이 주신 샘플 URL
    test_url = "https://www.creative-biolabs.com/bsab/bispecific-ang2-vegfa-tandem-diabody-12729.htm"
    
    logger.info(f"🚀 [TEST] Starting direct save test for: {test_url}")
    
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        context = await browser.new_context(user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36")
        page = await context.new_page()
        
        try:
            await page.goto(test_url, timeout=60000)
            await asyncio.sleep(5) # 충분한 로딩 대기
            
            title = await page.title()
            product_name = title.split('|')[0].strip()
            body_text = await page.inner_text("body")
            
            # 1. 상세 스펙 파싱 (Table)
            specs = await page.evaluate("""
                () => {
                    const data = {};
                    document.querySelectorAll('table tr').forEach(row => {
                        const cells = row.querySelectorAll('th, td');
                        if (cells.length >= 2) {
                            data[cells[0].innerText.trim().replace(/:$/, '')] = cells[1].innerText.trim();
                        }
                    });
                    return data;
                }
            """)
            
            # 2. 유전 정보 추출 (Regex)
            # Target 1Ang2 Gene ID285 UniProt IDO15123 ...
            # Target 2VEGFA Gene ID7422 UniProt IDP15692 ...
            
            def extract_target_details(text, index):
                # Target X 기호와 다음 Target 사이의 텍스트 추출
                pattern = f"Target {index}(.*?)Target {index+1}" if index == 1 else f"Target {index}(.*?)Alternative Names"
                match = re.search(pattern, text, re.DOTALL | re.IGNORECASE)
                if not match and index == 2: # 마지막 타겟일 경우
                     match = re.search(f"Target {index}(.*)", text, re.DOTALL | re.IGNORECASE)
                
                if match:
                    chunk = match.group(1)
                    symbol = chunk.split("Gene ID")[0].strip() if "Gene ID" in chunk else ""
                    gene_id = re.search(r"Gene ID\s*(\d+)", chunk, re.IGNORECASE)
                    uniprot_id = re.search(r"UniProt ID\s*([A-Z0-9]+)", chunk, re.IGNORECASE)
                    return {
                        "symbol": symbol,
                        "gene_id": gene_id.group(1) if gene_id else None,
                        "uniprot_id": uniprot_id.group(1) if uniprot_id else None
                    }
                return None

            t1 = extract_target_details(body_text, 1)
            t2 = extract_target_details(body_text, 2)
            targets = [t for t in [t1, t2] if t and t["symbol"]]
            
            logger.info(f"Parsed Targets: {targets}")

            # 3. DB 저장 시작
            # A. Antibody Library
            cat_no = f"CB-{test_url.split('-')[-1].replace('.htm', '')}"
            ab_data = {
                "product_name": product_name,
                "cat_no": cat_no,
                "antibody_format": specs.get("Type"),
                "host_species": specs.get("Host Animal 1"),
                "isotype": specs.get("Isotype"),
                "related_disease": specs.get("Related Disease"),
                "full_spec": specs,
                "source_url": test_url,
                "summary": f"Bispecific antibody targeting {specs.get('Targets')}"
            }
            
            logger.info("Saving to antibody_library...")
            ab_res = supabase.table("antibody_library").upsert(ab_data, on_conflict="cat_no").execute()
            if not ab_res.data:
                logger.error("❌ Failed to save antibody_library")
                return
            
            ab_id = ab_res.data[0]['id']
            logger.info(f"✅ Antibody Saved. ID: {ab_id}")

            # B. Target Master & Map
            for i, t in enumerate(targets):
                logger.info(f"Saving Target: {t['symbol']}...")
                t_master = {
                    "target_symbol": t["symbol"],
                    "gene_id": int(t["gene_id"]) if t["gene_id"] else None,
                    "uniprot_id": t["uniprot_id"]
                }
                t_res = supabase.table("target_master").upsert(t_master, on_conflict="target_symbol").execute()
                if t_res.data:
                    t_id = t_res.data[0]['id']
                    # Map
                    supabase.table("application_map").insert({
                        "antibody_id": ab_id,
                        "target_id": t_id,
                        "target_role": f"Target {i+1}",
                        "interaction_type": "Primary Binder"
                    }).execute()
                    logger.info(f"✅ Target {t['symbol']} mapped.")

            logger.info("🎉 [SUCCESS] Single product save test complete!")

        except Exception as e:
            logger.error(f"❌ Test Failed: {e}")
        finally:
            await browser.close()

if __name__ == "__main__":
    asyncio.run(test_single_product())
