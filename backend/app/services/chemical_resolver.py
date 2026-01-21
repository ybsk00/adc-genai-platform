import pubchempy as pcp
import logging
import time
from threading import Lock
from typing import Optional, Dict, Any

logger = logging.getLogger(__name__)

# RDKit 가용성 체크
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logger.warning("⚠️ RDKit not installed. SMILES validation will be skipped.")


class ChemicalResolver:
    """
    약물 이름을 기반으로 PubChem에서 SMILES 및 화학 정보를 조회하는 서비스
    Rate limiting + RDKit 검증으로 안전하게 처리
    """
    
    # Rate limiting 설정
    _last_call = 0
    _lock = Lock()
    MIN_INTERVAL = 1.0  # 최소 1초 간격
    
    @staticmethod
    def _wait_for_rate_limit():
        """Rate limit 대기"""
        with ChemicalResolver._lock:
            elapsed = time.time() - ChemicalResolver._last_call
            if elapsed < ChemicalResolver.MIN_INTERVAL:
                wait_time = ChemicalResolver.MIN_INTERVAL - elapsed
                logger.debug(f"⏳ PubChem rate limit: waiting {wait_time:.2f}s")
                time.sleep(wait_time)
            ChemicalResolver._last_call = time.time()
    
    @staticmethod
    def fetch_safe_smiles(drug_name: str) -> Dict[str, Any]:
        """
        안전하게 PubChem에서 SMILES를 가져오는 방탄 함수
        RDKit으로 화학적 무결성 검증
        """
        result = {"smiles": None, "status": "NOT_FOUND", "mw": 0, "cid": None}
        
        if not drug_name or drug_name.lower() in ["unknown", "n/a", "none"]:
            return result
        
        # Rate limit 적용
        ChemicalResolver._wait_for_rate_limit()
        
        try:
            # 1. PubChem 검색
            compounds = pcp.get_compounds(drug_name, 'name')
            
            if not compounds:
                logger.warning(f"⚠️ PubChem: '{drug_name}' 검색 결과 없음.")
                return result
            
            # 2. 첫 번째 후보 가져오기
            candidate = compounds[0]
            candidate_smiles = candidate.isomeric_smiles
            result["cid"] = candidate.cid
            
            # 3. RDKit 검증 (가능한 경우)
            if RDKIT_AVAILABLE and candidate_smiles:
                mol = Chem.MolFromSmiles(candidate_smiles)
                if mol:
                    mw = Descriptors.MolWt(mol)
                    result["mw"] = round(mw, 2)
                    
                    # ADC Payload용 분자량 필터 (너무 작거나 크면 의심)
                    if 200 < mw < 2000:
                        result["smiles"] = candidate_smiles
                        result["status"] = "VERIFIED"
                        logger.info(f"✅ SMILES verified for {drug_name}: MW={mw:.1f}")
                    else:
                        result["status"] = "REVIEW_NEEDED"
                        result["smiles"] = candidate_smiles
                        logger.warning(f"⚠️ MW out of range for {drug_name}: {mw:.1f}")
                else:
                    result["status"] = "INVALID_STRUCTURE"
                    logger.warning(f"❌ Invalid SMILES structure for {drug_name}")
            else:
                # RDKit 없으면 검증 없이 저장
                result["smiles"] = candidate_smiles
                result["status"] = "UNVERIFIED"
                logger.info(f"📝 SMILES fetched (unverified) for {drug_name}")
                
        except Exception as e:
            logger.error(f"❌ PubChem Error ({drug_name}): {str(e)}")
            result["status"] = "API_ERROR"
        
        return result
    
    @staticmethod
    def fetch_verified_smiles(drug_name: str) -> Optional[str]:
        """
        기존 호환성을 위한 래퍼 함수
        SMILES 문자열만 반환 (None 또는 문자열)
        """
        result = ChemicalResolver.fetch_safe_smiles(drug_name)
        return result["smiles"]

    @staticmethod
    def get_cid(drug_name: str) -> Optional[int]:
        """PubChem CID 조회"""
        result = ChemicalResolver.fetch_safe_smiles(drug_name)
        return result["cid"]


chemical_resolver = ChemicalResolver()

