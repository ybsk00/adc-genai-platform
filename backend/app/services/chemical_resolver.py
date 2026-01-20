import pubchempy as pcp
import logging
from typing import Optional

logger = logging.getLogger(__name__)

class ChemicalResolver:
    """
    약물 이름을 기반으로 PubChem에서 SMILES 및 화학 정보를 조회하는 서비스
    """
    
    @staticmethod
    def fetch_verified_smiles(drug_name: str) -> Optional[str]:
        """
        PubChem에서 정확한 SMILES 코드를 찾아옵니다.
        """
        if not drug_name or drug_name.lower() == "unknown":
            return None
            
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
        try:
            compounds = pcp.get_compounds(drug_name, 'name')
            return compounds[0].cid if compounds else None
        except:
            return None

chemical_resolver = ChemicalResolver()
