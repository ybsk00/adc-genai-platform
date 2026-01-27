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
    # Non-critical exit, just won't validate MW
    pass

# Load env
load_dotenv()

# Setup Logging
logger = logging.getLogger("AI_Linker_Refiner")
logger.setLevel(logging.INFO)
if logger.hasHandlers(): logger.handlers.clear()
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(ch)

# Supabase & Gemini
url: str = os.environ.get("SUPABASE_URL")
key: str = os.environ.get("SUPABASE_SERVICE_KEY")
supabase: Client = create_client(url, key)

genai.configure(api_key=os.getenv("GOOGLE_API_KEY"))
model = genai.GenerativeModel('gemini-2.0-flash') # Using newer model for better chemistry understanding

# Known Linkers Dictionary (Can be expanded)
KNOWN_LINKERS = {
    "SMCC": "C1CCC(CC1)C(=O)N2CCC(=O)N(C2=O)CC3=CC=C(C=C3)N4C(=O)C=CC4=O",
    "mc-Val-Cit-PAB": "CC(C)[C@@H](NC(=O)CCCCC1=CC(=O)N(C1)CC(=O)O)C(=O)N[C@@H](CCCNC(=O)N)C(=O)NCC2=CC=C(C=C2)CO",
    "vc-PAB": "CC(C)[C@@H](N)C(=O)N[C@@H](CCCNC(=O)N)C(=O)NCC1=CC=C(C=C1)CO", # Val-Cit-PAB
    "MMAE Linker": "CC(C)[C@@H](NC(=O)CCCCC1=CC(=O)N(C1)CC(=O)O)C(=O)N[C@@H](CCCNC(=O)N)C(=O)NCC2=CC=C(C=C2)COC(=O)N(C)C", # mc-Val-Cit-PABC
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

async def ask_llm_for_cid(name: str, cas: str, summary: str):
    """LLM에게 SMILES 대신 CID를 물어봄 (검색 유도)"""
    clean_name = name.split('(')[0].split('CAT#:')[0].strip()
    
    prompt = f"""
You are an expert in Antibody-Drug Conjugate (ADC) chemistry. 
Find the PubChem CID (Compound ID) for the following ADC Linker.

    Product Name: {clean_name}
    CAS: {cas}
    Description: {summary}

    Instruction:
    - Search for the PubChem CID.
    - If it's a common linker like SMCC, SPDB, mc-Val-Cit-PAB, identify the CID for the linker molecule itself.
    - Do NOT generate a SMILES code directly unless you are 100% sure. CID is preferred. 
    
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
    """RDKit으로 SMILES 유효성 및 분자량 오차 검증"""
    if 'rdkit' not in globals():
        return True, "RDKit skipped", smiles

    if not smiles or smiles == "None":
        return False, "SMILES is None or empty", None

    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if not mol:
            return False, "Invalid SMILES syntax (RDKit parsing failed)", None
        
        calc_mw = Descriptors.MolWt(mol)
        
        if not target_mw_str:
            return True, f"Valid Structure (Target MW missing, Calc: {calc_mw:.2f})", smiles

        # Clean target MW
        target_mw_clean = re.sub(r'[^0-9.]', '', str(target_mw_str))
        if not target_mw_clean:
            return True, f"Valid Structure (Target MW unreadable, Calc: {calc_mw:.2f})", smiles
            
        target_mw = float(target_mw_clean)
        
        # 1차 비교: 원본 그대로 비교
        error = abs(calc_mw - target_mw)
        if error < 10.0: # Linkers might have salts or variations, allow 10.0
            return True, f"MW Match (Calc: {calc_mw:.2f}, Target: {target_mw}, Error: {error:.2f})", smiles

        # 2차 비교 (패자부활전): Desalt
        remover = SaltRemover()
        mol_desalted = remover.StripMol(mol, dontRemoveEverything=True)
        calc_mw_desalted = Descriptors.MolWt(mol_desalted)
        
        error_desalted = abs(calc_mw_desalted - target_mw)
        
        if error_desalted < 10.0:
            smiles_desalted = Chem.MolToSmiles(mol_desalted)
            return True, f"Desalted Match (Orig: {calc_mw:.2f} -> Desalt: {calc_mw_desalted:.2f}, Target: {target_mw})", smiles_desalted
        
        return False, f"MW Mismatch (Calc: {calc_mw:.2f}/Desalt: {calc_mw_desalted:.2f} vs Target: {target_mw})", None
            
    except Exception as e:
        return False, f"RDKit Error: {e}", None

async def run_linker_refinement():
    logger.info("🔗 Starting Linker Data Refinement Pipeline...")
    
    # 1. Fetch unrefined Linkers
    res = supabase.table("commercial_reagents")\
        .select("*")\
        .eq("source_name", "Creative Biolabs")\
        .eq("category", "ADC Linkers")\
        .is_("smiles_code", "null")\
        .execute()
    
    targets = res.data
    logger.info(f"🎯 Linker Targets found: {len(targets)}")
    
    for item in targets:
        rid = item['id']
        name = item['product_name']
        cas = item['cas_number']
        target_mw = item.get('molecular_weight')
        summary = item.get('summary', '') or ''
        
        logger.info(f"🔬 Analyzing: {name} (CAS: {cas}, MW: {target_mw})")
        
        candidate = None
        
        # Strategy 0: Known Linkers (Dictionary)
        for key, gold_smiles in KNOWN_LINKERS.items():
            if key.lower() in name.lower():
                candidate = {"smiles": gold_smiles, "source": "Internal Dictionary", "mw": None}
                logger.info(f"   💡 Found in Internal Dictionary: {key}")
                break

        # Strategy 1: PubChem by CAS
        if not candidate and cas:
            candidate = await fetch_pubchem(cas, 'name')
        
        # Strategy 2: PubChem by Name
        if not candidate:
            candidate = await fetch_pubchem(name, 'name')
            
        # Strategy 3: LLM CID Search
        if not candidate:
            logger.info("   ⚠️ Direct API failed. Asking Gemini for CID...")
            cid = await ask_llm_for_cid(name, cas, summary)
            if cid:
                logger.info(f"   🤖 Gemini found CID: {cid}. Fetching from PubChem...")
                candidate = await fetch_pubchem_by_cid(cid)
        
        # Validation & Update
        if candidate:
            smiles = candidate['smiles']
            is_valid, reason, final_smiles = validate_rdkit(smiles, target_mw)
            
            if candidate['source'] == "Internal Dictionary":
                is_valid = True # Trust internal dictionary
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
                    "summary": f"Refined via {candidate['source']} | {summary}"[:500] # Append original summary
                }).eq("id", rid).execute()
            else:
                logger.warning(f"   ❌ REJECTED: {reason}")
                # Log the failure so we don't retry endlessly (optional, or just leave as null to retry later)
                supabase.table("commercial_reagents").update({
                    "summary": f"Refinement Failed: {reason} | {summary}"[:500]
                }).eq("id", rid).execute()
        else:
            logger.warning("   🚫 No structure found.")
            
        await asyncio.sleep(1)

if __name__ == "__main__":
    asyncio.run(run_linker_refinement())
