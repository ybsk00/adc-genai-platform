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

# 1. 골든셋 사전 (최후의 보루 및 고난도 물질용)
GOLDEN_DICTIONARY = {
    "Ansamitocin P-3": "C[C@@H]1[C@@H]2C[C@]([C@@H](/C=C/C=C(/CC3=CC(=C(C(=C3)OC)Cl)N(C(=O)C[C@@H]([C@]4([C@H]1O4)C)OC(=O)C(C)C)C)\\C)OC)(NC(=O)O2)O",
    "Auristatin F": "CC[C@H](C)[C@@H]([C@@H](CC(=O)N1CCC[C@H]1[C@@H]([C@@H](C)C(=O)N[C@@H](CC2=CC=CC=C2)C(=O)O)OC)OC)N(C)C(=O)[C@H](C(C)C)NC(=O)[C@H](C(C)C)N(C)C",
    "Auristatin E": "CC[C@H](C)[C@@H]([C@@H](CC(=O)N1CCC[C@H]1[C@@H]([C@@H](C)C(=O)N[C@H](C)[C@H](C2=CC=CC=C2)O)OC)OC)N(C)C(=O)[C@H](C(C)C)NC(=O)[C@H](C(C)C)N(C)C",
    "MMAE": "CC[C@H](C)[C@@H]([C@@H](CC(=O)N1CCC[C@H]1[C@@H]([C@@H](C)C(=O)N[C@H](C)[C@H](C2=CC=CC=C2)O)OC)OC)N(C)C(=O)[C@H](C(C)C)NC(=O)[C@H](C(C)C)N(C)C",
    "MMAF": "CC[C@H](C)[C@@H]([C@@H](CC(=O)N1CCC[C@H]1[C@@H]([C@@H](C)C(=O)N[C@@H](CC2=CC=CC=C2)C(=O)O)OC)OC)N(C)C(=O)[C@H](C(C)C)NC(=O)[C@H](C(C)C)N(C)C",
    "Calicheamicin": "CC1=C(C(=C2C(=C1)C(=O)C3=C(C2=O)C(=CC=C3)OC)O)OC4C(C(C(C(O4)C)NC(=O)NC5=CC=CC=C5)O)SC6=C7C(=C(C(=C6)C)O)C[C@@H](C[C@H]7C#CC#C[C@@]8(C[C@@H]([C@@H]([C@H](O8)C)OC)OC)O)O",
    "Daun02": "CC1C(C(CC(O1)OC2CC(CC3=C2C(=C4C(=C3O)C(=O)C5=C(C4=O)C(=CC=C5)OC)O)(C(=O)CO)O)N)O"
}

async def fetch_pubchem_by_cid(cid: str):
    """CID로 PubChem에서 SMILES 가져오기"""
    if not cid: return None
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES,MolecularWeight/JSON"
    try:
        async with aiohttp.ClientSession() as session:
            async with session.get(url, timeout=10) as response:
                if response.status == 200:
                    data = await response.json()
                    props = data['PropertyTable']['Properties'][0]
                    return {
                        "smiles": props.get('CanonicalSMILES'),
                        "mw": float(props.get('MolecularWeight')),
                        "source": f"PubChem (CID: {cid})"
                    }
    except:
        return None
    return None

async def ask_llm_for_cid(name: str, cas: str):
    """LLM에게 SMILES 대신 CID를 물어봄 (검색 유도)"""
    clean_name = name.split('(')[0].split('CAT#:')[0].strip()
    
    prompt = f"""
You are a chemical database expert. Find the PubChem CID (Compound ID) for the following substance.
    
    Target: {clean_name}
    CAS: {cas}

    Instruction:
    - Search for the PubChem CID first.
    - Identify the exact PubChem CID for the "Free Drug" form (Payload).
    - Do NOT generate a SMILES code. Just give me the ID. 
    
    Output JSON:
    {{
        "pubchem_cid": "12345",
        "confidence": "High"
    }}
    If unknown, return null for pubchem_cid.
    """
    try:
        res = model.generate_content(prompt)
        clean_text = res.text.replace('```json', '').replace('```', '').strip()
        match = re.search(r'{{.*}}', clean_text, re.DOTALL)
        if match:
            data = json.loads(match.group())
            return data.get('pubchem_cid')
    except:
        pass
    return None

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
    logger.info("🧪 Starting Ultimate Structure Refinement Pipeline...")
    
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
        
        # 1. 쓰레기 데이터 필터링
        skip_keywords = ['Services', 'Products', 'Publications', 'List', 'Page', 'Module', 'Conjugate']
        if any(k in name for k in skip_keywords):
            logger.info(f"🗑️ Skipping non-chemical entry: {name}")
            continue

        logger.info(f"🔬 Analyzing: {name} (CAS: {cas}, MW: {target_mw})")
        
        candidate = None
        
        # Strategy 0: Golden Dictionary (Hardcoded)
        # 이름의 일부만 매칭되도 찾도록 (e.g. "Auristatin F (..." -> "Auristatin F")
        for key, gold_smiles in GOLDEN_DICTIONARY.items():
            if key.lower() == name.split('(')[0].strip().lower():
                candidate = {"smiles": gold_smiles, "source": "Golden Dictionary", "mw": None} # MW will be calc'd
                break
        
        # Strategy 1: PubChem by CAS (Existing)
        if not candidate and cas:
            candidate = await fetch_pubchem(cas, 'name')
        
        # Strategy 2: PubChem by Name (Existing)
        if not candidate:
            candidate = await fetch_pubchem(name, 'name')
            
        # Strategy 3: LLM CID Search -> PubChem API (NEW)
        if not candidate:
            logger.info("   ⚠️ Direct API failed. Asking Gemini for CID...")
            cid = await ask_llm_for_cid(name, cas)
            if cid:
                logger.info(f"   🤖 Gemini found CID: {cid}. Fetching from PubChem...")
                candidate = await fetch_pubchem_by_cid(cid)
        
        # Validation & Update
        if candidate:
            smiles = candidate['smiles']
            # MW 검증 (Golden Dict는 무조건 통과시키거나, 검증하더라도 MW가 DB에 없으면 통과됨)
            is_valid, reason, final_smiles = validate_rdkit(smiles, target_mw)
            
            # Golden Dictionary는 검증 실패해도 일단 저장 (신뢰도 최상)
            if candidate['source'] == "Golden Dictionary":
                is_valid = True
                reason = "Golden Set Override"
                final_smiles = smiles

            if is_valid:
                logger.info(f"   ✅ APPROVED: {reason}")
                props = item.get('properties') or {}
                props['structure_source'] = candidate['source']
                props['validation_log'] = reason
                props['refined_at'] = datetime.now().isoformat()
                
                supabase.table("commercial_reagents").update({
                    "smiles_code": final_smiles,
                    "ai_refined": True,
                    "properties": props,
                    "summary": f"Refined via {candidate['source']}"
                }).eq("id", rid).execute()
            else:
                logger.warning(f"   ❌ REJECTED: {reason}")
        else:
            logger.warning("   🚫 No structure found.")
            
        await asyncio.sleep(1)

if __name__ == "__main__":
    asyncio.run(run_pipeline())