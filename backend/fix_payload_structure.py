import os
import asyncio
import logging
import json
import re
import aiohttp
from datetime import datetime
from dotenv import load_dotenv
import google.generativeai as genai
from supabase import create_client, Client

# Load env
load_dotenv()

# Setup Logging
logger = logging.getLogger("AI_Structure_Fixer")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(ch)

# Supabase Setup
url = os.getenv("SUPABASE_URL")
key = os.getenv("SUPABASE_SERVICE_KEY")
supabase: Client = create_client(url, key)

# Gemini Setup
genai.configure(api_key=os.getenv("GOOGLE_API_KEY"))
model = genai.GenerativeModel('gemini-1.5-flash')

async def fetch_smiles_from_pubchem(cas_number: str):
    """CAS 번호로 PubChem에서 SMILES 가져오기"""
    if not cas_number or cas_number.lower() == 'none':
        return None
    
    # Clean CAS number (sometimes has spaces or prefix)
    cas_clean = re.sub(r'[^0-9-]', '', cas_number)
    
    pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas_clean}/property/CanonicalSMILES,MolecularFormula,MolecularWeight/JSON"
    
    try:
        async with aiohttp.ClientSession() as session:
            async with session.get(pubchem_url, timeout=10) as response:
                if response.status == 200:
                    data = await response.json()
                    props = data['PropertyTable']['Properties'][0]
                    return {
                        "smiles": props.get('CanonicalSMILES'),
                        "formula": props.get('MolecularFormula'),
                        "mw": str(props.get('MolecularWeight'))
                    }
                else:
                    # Try by Name if CAS fails (Optional)
                    return None
    except Exception as e:
        logger.error(f"PubChem API Error ({cas_number}): {e}")
        return None

async def verify_with_ai(product_name: str, smiles: str, formula: str):
    """AI를 사용하여 제품명과 구조식의 일치 여부 검증"""
    prompt = f"""
    Check if the following chemical information is consistent.
    Product Name: {product_name}
    Proposed SMILES: {smiles}
    Proposed Formula: {formula}

    Does this SMILES/Formula correctly represent the product name? 
    Answer in JSON format:
    {{
        "is_match": true/false,
        "reason": "short explanation",
        "confidence": 0.0 to 1.0
    }}
    """
    
    try:
        response = model.generate_content(prompt)
        # Extract JSON from response
        match = re.search(r'{{.*}}', response.text, re.DOTALL)
        if match:
            return json.loads(match.group())
        return None
    except Exception as e:
        logger.error(f"AI Verification Error: {e}")
        return None

async def run_fixer():
    logger.info("🛠️ Starting AI Structure Fixer...")
    
    # 1. 대상 레코드 조회 (CAS는 있는데 SMILES는 없는 것들)
    res = supabase.table("commercial_reagents")\
        .select("*")\
        .not_.is_("cas_number", "null")\
        .is_("smiles_code", "null")\
        .eq("source_name", "Creative Biolabs")\
        .execute()
    
    targets = res.data
    logger.info(f"🔎 Found {len(targets)} records to refine.")
    
    for item in targets:
        record_id = item['id']
        name = item['product_name']
        cas = item['cas_number']
        
        logger.info(f"🔄 Processing: {name} (CAS: {cas})")
        
        # 2. PubChem에서 정보 가져오기
        pc_data = await fetch_smiles_from_pubchem(cas)
        
        if pc_data and pc_data['smiles']:
            smiles = pc_data['smiles']
            formula = pc_data['formula']
            
            # 3. AI 검증
            v_res = await verify_with_ai(name, smiles, formula)
            
            if v_res and v_res.get('is_match'):
                logger.info(f"✅ AI Verified: {name} matches {smiles}")
                
                # 4. DB 업데이트
                supabase.table("commercial_reagents").update({
                    "smiles_code": smiles,
                    "formula": formula,
                    "ai_refined": True,
                    "summary": f"Structure verified via PubChem and AI. (Reason: {v_res.get('reason')})"
                }).eq("id", record_id).execute()
            else:
                reason = v_res.get('reason') if v_res else "Verification failed"
                logger.warning(f"❌ AI Rejected: {name} (Reason: {reason})")
                supabase.table("commercial_reagents").update({
                    "ai_refined": False,
                    "summary": f"AI Refinement Failed: {reason}"
                }).eq("id", record_id).execute()
        else:
            logger.warning(f"⚠️ No data found in PubChem for CAS: {cas}")
        
        # Rate limit safety
        await asyncio.sleep(2)

if __name__ == "__main__":
    asyncio.run(run_fixer())
