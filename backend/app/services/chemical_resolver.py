import pubchempy as pcp
import logging
import time
from threading import Lock
from typing import Optional

logger = logging.getLogger(__name__)

class ChemicalResolver:
    """
    약물 이름을 기반으로 PubChem에서 SMILES 및 화학 정보를 조회하는 서비스
    Rate limiting으로 API 차단 방지
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
    def fetch_verified_smiles(drug_name: str) -> Optional[str]:
        """
        PubChem에서 정확한 SMILES 코드를 찾아옵니다.
        """
        if not drug_name or drug_name.lower() == "unknown":
            return None
        
        # Rate limit 적용
        ChemicalResolver._wait_for_rate_limit()
            
        try:
            # 1. 이름으로 화합물 검색
            compounds = pcp.get_compounds(drug_name, 'name')
            
            if compounds:
                # 2. 가장 관련성 높은 첫 번째 결과의 Isomeric SMILES 반환
                smiles = compounds[0].isomeric_smiles
                logger.info(f"✅ PubChem found SMILES for {drug_name}: {smiles}")
                return smiles
            else:
                logger.warning(f"❓ No PubChem entry found for {drug_name}")
                return None
                
        except Exception as e:
            logger.error(f"❌ PubChem API Error for {drug_name}: {e}")
            return None

    @staticmethod
    def get_cid(drug_name: str) -> Optional[int]:
        """PubChem CID 조회"""
        ChemicalResolver._wait_for_rate_limit()
        try:
            compounds = pcp.get_compounds(drug_name, 'name')
            return compounds[0].cid if compounds else None
        except:
            return None

chemical_resolver = ChemicalResolver()
