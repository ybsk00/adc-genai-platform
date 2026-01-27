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

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem.SaltRemover import SaltRemover
except ImportError:
    print("❌ RDKit is missing. Please install: pip install rdkit")
    exit(1)

# Load env
load_dotenv()

# Setup Logging
logger = logging.getLogger("AI_Structure_Fixer_Adv")
logger.setLevel(logging.INFO)
if logger.hasHandlers(): logger.handlers.clear()
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(ch)

# Supabase & Gemini
supabase: Client = create_client(os.getenv("SUPABASE_URL"), os.getenv("SUPABASE_SERVICE_KEY"))
genai.configure(api_key=os.getenv("GOOGLE_API_KEY"))
model = genai.GenerativeModel('gemini-2.0-flash')

async def fetch_pubchem(identifier: str, namespace: str = 'name'):
    """PubChem API Call helper"""
    if not identifier: return None
    clean_id = identifier.strip()
    if namespace == 'name':
        # 이름 세척
        clean_id = clean_id.split('(')[0].split('CAT#:')[0].strip()
        clean_id = re.sub(r'\s+', '%20', clean_id)
    
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{namespace}/{clean_id}/property/CanonicalSMILES,MolecularFormula,MolecularWeight/JSON"
    
    try:
        async with aiohttp.ClientSession() as session:
            async with session.get(url, timeout=10) as response:
                if response.status == 200:
                    data = await response.json()
                    if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                        props = data['PropertyTable']['Properties'][0]
                        
                        smiles = props.get('CanonicalSMILES')
                        if not smiles: return None
                        
                        mw_raw = props.get('MolecularWeight')
                        mw = float(mw_raw) if mw_raw else None
                        
                        return {
                            "smiles": smiles,
                            "mw": mw,
                            "source": f"PubChem ({namespace})"
                        }
    except Exception as e:
        return None
    return None

async def generate_smiles_via_llm(name: str, cas: str):
    """LLM에게 SMILES 생성을 요청"""
    # 이름 세척 (API와 동일하게 적용)
    clean_name = name.split('(')[0].split('CAT#:')[0].strip()
    
    prompt = f"""
    You are a chemical expert. Provide the Canonical SMILES for the following substance.
    Product Name: {clean_name}
    CAS Number: {cas}

    Output strictly in JSON format:
    {{
        "smiles": "INSERT_SMILES_HERE",
        "source": "Gemini 2.0 Flash Knowledge"
    }}
    If unknown, return null for smiles.
    """
    try:
        res = model.generate_content(prompt)
        match = re.search(r'{{.*}}', res.text, re.DOTALL)
        if match:
            data = json.loads(match.group())
            if data.get('smiles') and data['smiles'] != "null":
                return data
            else:
                # 디버그: 왜 null을 줬는지 확인
                logger.warning(f"   🤖 Gemini returned null/empty for '{clean_name}'. Raw: {data}")
        else:
             logger.warning(f"   🤖 Gemini response parse failed for '{clean_name}'. Raw: {res.text[:100]}...")
             
    except Exception as e:
        logger.error(f"LLM Error: {e}")
    return None

def validate_rdkit(smiles: str, target_mw_str: str):
    """RDKit으로 SMILES 유효성 및 분자량 오차 검증 (with Desalt Fallback)"""
    if not smiles or smiles == "None":
        return False, "SMILES is None or empty", None

    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if not mol:
            return False, "Invalid SMILES syntax (RDKit parsing failed)", None
        
        calc_mw = Descriptors.MolWt(mol)
        
        # 타겟 MW가 없으면 구조적 유효성만 체크
        if not target_mw_str:
            return True, f"Valid Structure (Target MW missing, Calc: {calc_mw:.2f})", smiles

        # Clean target MW
        target_mw_clean = re.sub(r'[^0-9.]', '', str(target_mw_str))
        if not target_mw_clean:
            return True, f"Valid Structure (Target MW unreadable, Calc: {calc_mw:.2f})", smiles
            
        target_mw = float(target_mw_clean)
        
        # 1차 비교: 원본 그대로 비교
        error = abs(calc_mw - target_mw)
        if error < 5.0: # 허용 오차 5.0
            return True, f"MW Match (Calc: {calc_mw:.2f}, Target: {target_mw}, Error: {error:.2f})", smiles

        # 2차 비교 (패자부활전): Desalt 후 비교
        remover = SaltRemover()
        mol_desalted = remover.StripMol(mol, dontRemoveEverything=True)
        calc_mw_desalted = Descriptors.MolWt(mol_desalted)
        
        error_desalted = abs(calc_mw_desalted - target_mw)
        
        if error_desalted < 5.0:
            smiles_desalted = Chem.MolToSmiles(mol_desalted)
            return True, f"Desalted Match (Orig: {calc_mw:.2f} -> Desalt: {calc_mw_desalted:.2f}, Target: {target_mw})", smiles_desalted
        
        # 최종 실패
        return False, f"MW Mismatch (Calc: {calc_mw:.2f}/Desalt: {calc_mw_desalted:.2f} vs Target: {target_mw})", None
            
    except Exception as e:
        return False, f"RDKit Error: {e}", None

async def run_pipeline():
    logger.info("🧪 Starting Advanced Structure Refinement Pipeline (with Desalter)...")
    
    # 1. 대상 조회
    res = supabase.table("commercial_reagents")\
        .select("*")\
        .eq("source_name", "Creative Biolabs")\
        .is_("smiles_code", "null")\
        .execute()
    
    targets = res.data
    logger.info(f"🎯 Targets found: {len(targets)}")
    
    for item in targets:
        rid = item['id']
        name = item['product_name']
        cas = item['cas_number']
        target_mw = item.get('molecular_weight')
        
        logger.info(f"🔬 Analyzing: {name} (CAS: {cas}, MW: {target_mw})")
        
        candidate = None
        
        # Step 1: PubChem by CAS
        if cas:
            candidate = await fetch_pubchem(cas, 'name')
        
        # Step 2: PubChem by Name
        if not candidate:
            candidate = await fetch_pubchem(name, 'name')
            
        # Step 3: LLM Fallback
        if not candidate:
            logger.info("   ⚠️ API failed. Asking Gemini...")
            llm_res = await generate_smiles_via_llm(name, cas)
            if llm_res:
                candidate = {
                    "smiles": llm_res['smiles'],
                    "source": "Gemini 2.0 Flash Knowledge"
                }
        
        # Validation & Update
        if candidate:
            smiles = candidate['smiles']
            # validate_rdkit returns (is_valid, reason, final_smiles)
            is_valid, reason, final_smiles = validate_rdkit(smiles, target_mw)
            
            if is_valid:
                logger.info(f"   ✅ APPROVED: {reason}")
                
                # Update DB
                props = item.get('properties') or {}
                props['structure_source'] = candidate['source']
                props['validation_log'] = reason
                props['refined_at'] = datetime.now().isoformat()
                
                supabase.table("commercial_reagents").update({
                    "smiles_code": final_smiles,
                    "ai_refined": True,
                    "properties": props,
                    "summary": f"Structure refined via {candidate['source']}. {reason}"
                }).eq("id", rid).execute()
            else:
                logger.warning(f"   ❌ REJECTED: {reason}")
        else:
            logger.warning("   🚫 No structure found from any source.")
            
        await asyncio.sleep(1)

if __name__ == "__main__":
    asyncio.run(run_pipeline())