"""
PubChem Service - CAS/Name to SMILES Lookup with MW Validation
Features:
- CAS Number to SMILES lookup via PubChem API
- Desalt + MW Validation using RDKit
- In-memory caching to prevent API overload
- Rate limiting protection
"""
import asyncio
import httpx
from typing import Optional, Dict, Any, Tuple
from datetime import datetime, timedelta
from functools import lru_cache
import logging

logger = logging.getLogger("PubChem_Service")

# Simple in-memory cache with TTL
class SimpleCache:
    def __init__(self, ttl_minutes: int = 60):
        self._cache: Dict[str, Tuple[Any, datetime]] = {}
        self._ttl = timedelta(minutes=ttl_minutes)

    def get(self, key: str) -> Optional[Any]:
        if key in self._cache:
            value, timestamp = self._cache[key]
            if datetime.now() - timestamp < self._ttl:
                return value
            else:
                del self._cache[key]
        return None

    def set(self, key: str, value: Any):
        self._cache[key] = (value, datetime.now())

    def clear(self):
        self._cache.clear()


class PubChemService:
    """
    PubChem API Service with caching and MW validation
    """
    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    def __init__(self):
        self._cache = SimpleCache(ttl_minutes=120)  # 2시간 캐시
        self._rate_limit_delay = 0.3  # 300ms between requests (PubChem guideline)
        self._last_request_time = datetime.min

    async def _rate_limit(self):
        """Rate limiting to avoid PubChem API overload"""
        elapsed = (datetime.now() - self._last_request_time).total_seconds()
        if elapsed < self._rate_limit_delay:
            await asyncio.sleep(self._rate_limit_delay - elapsed)
        self._last_request_time = datetime.now()

    async def get_smiles_by_cas(self, cas_number: str) -> Dict[str, Any]:
        """
        CAS Number로 SMILES 조회
        Returns: {smiles_code, molecular_weight, formula, pubchem_cid, validation_status, error}
        """
        if not cas_number:
            return {"error": "CAS number is required"}

        # Normalize CAS
        cas_clean = cas_number.strip().replace(" ", "")

        # Check cache
        cache_key = f"cas:{cas_clean}"
        cached = self._cache.get(cache_key)
        if cached:
            logger.info(f"Cache hit for CAS: {cas_clean}")
            return cached

        try:
            await self._rate_limit()

            async with httpx.AsyncClient(timeout=30.0) as client:
                # Step 1: CAS -> CID
                cid_url = f"{self.BASE_URL}/compound/name/{cas_clean}/cids/JSON"
                cid_res = await client.get(cid_url)

                if cid_res.status_code != 200:
                    return {"error": f"CAS not found in PubChem: {cas_clean}"}

                cid_data = cid_res.json()
                cid = cid_data.get("IdentifierList", {}).get("CID", [None])[0]

                if not cid:
                    return {"error": "No CID found for this CAS"}

                # Step 2: CID -> Properties (SMILES, MW, Formula)
                await self._rate_limit()
                props_url = f"{self.BASE_URL}/compound/cid/{cid}/property/CanonicalSMILES,MolecularWeight,MolecularFormula/JSON"
                props_res = await client.get(props_url)

                if props_res.status_code != 200:
                    return {"error": "Failed to fetch compound properties"}

                props_data = props_res.json()
                props = props_data.get("PropertyTable", {}).get("Properties", [{}])[0]

                smiles = props.get("CanonicalSMILES")
                mw = props.get("MolecularWeight")
                formula = props.get("MolecularFormula")

                if not smiles:
                    return {"error": "No SMILES found for this compound"}

                # Step 3: Desalt + MW Validation (RDKit)
                validation_result = await self._validate_smiles_with_rdkit(smiles, mw)

                result = {
                    "smiles_code": smiles,
                    "molecular_weight": mw,
                    "formula": formula,
                    "pubchem_cid": cid,
                    "validation_status": validation_result.get("status"),
                    "desalted_smiles": validation_result.get("desalted_smiles"),
                    "calculated_mw": validation_result.get("calculated_mw"),
                    "mw_difference": validation_result.get("mw_difference"),
                    "source": "PubChem"
                }

                # Cache result
                self._cache.set(cache_key, result)
                logger.info(f"Fetched and cached SMILES for CAS: {cas_clean}")

                return result

        except httpx.TimeoutException:
            return {"error": "PubChem API timeout"}
        except Exception as e:
            logger.error(f"PubChem API error: {e}")
            return {"error": f"API error: {str(e)}"}

    async def get_smiles_by_name(self, compound_name: str) -> Dict[str, Any]:
        """
        화합물명으로 SMILES 조회
        """
        if not compound_name:
            return {"error": "Compound name is required"}

        name_clean = compound_name.strip()

        # Check cache
        cache_key = f"name:{name_clean.lower()}"
        cached = self._cache.get(cache_key)
        if cached:
            return cached

        try:
            await self._rate_limit()

            async with httpx.AsyncClient(timeout=30.0) as client:
                # Name -> CID -> Properties
                cid_url = f"{self.BASE_URL}/compound/name/{name_clean}/cids/JSON"
                cid_res = await client.get(cid_url)

                if cid_res.status_code != 200:
                    return {"error": f"Compound not found: {name_clean}"}

                cid_data = cid_res.json()
                cid = cid_data.get("IdentifierList", {}).get("CID", [None])[0]

                if not cid:
                    return {"error": "No CID found"}

                await self._rate_limit()
                props_url = f"{self.BASE_URL}/compound/cid/{cid}/property/CanonicalSMILES,MolecularWeight,MolecularFormula/JSON"
                props_res = await client.get(props_url)

                props_data = props_res.json()
                props = props_data.get("PropertyTable", {}).get("Properties", [{}])[0]

                smiles = props.get("CanonicalSMILES")
                mw = props.get("MolecularWeight")
                formula = props.get("MolecularFormula")

                validation_result = await self._validate_smiles_with_rdkit(smiles, mw)

                result = {
                    "smiles_code": smiles,
                    "molecular_weight": mw,
                    "formula": formula,
                    "pubchem_cid": cid,
                    "validation_status": validation_result.get("status"),
                    "desalted_smiles": validation_result.get("desalted_smiles"),
                    "calculated_mw": validation_result.get("calculated_mw"),
                    "mw_difference": validation_result.get("mw_difference"),
                    "source": "PubChem"
                }

                self._cache.set(cache_key, result)
                return result

        except Exception as e:
            return {"error": f"API error: {str(e)}"}

    async def _validate_smiles_with_rdkit(self, smiles: str, target_mw: Optional[float] = None) -> Dict[str, Any]:
        """
        RDKit을 사용한 SMILES 검증 + Desalt + MW 비교
        - Desalt: 염(salt) 제거 후 주요 분자만 추출
        - MW 검증: Target MW와 ±1.0 오차 이내인지 확인
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors
            from rdkit.Chem.SaltRemover import SaltRemover

            # Parse SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {"status": "invalid", "error": "Invalid SMILES structure"}

            # Desalt (염 제거)
            remover = SaltRemover()
            desalted_mol = remover.StripMol(mol)

            # 가장 큰 fragment 선택 (염이 완전히 제거 안 됐을 경우)
            frags = Chem.GetMolFrags(desalted_mol, asMols=True)
            if frags:
                # 가장 무거운 fragment 선택
                desalted_mol = max(frags, key=lambda x: Descriptors.MolWt(x))

            desalted_smiles = Chem.MolToSmiles(desalted_mol)
            calculated_mw = round(Descriptors.MolWt(desalted_mol), 2)

            # MW 비교 (±1.0 허용)
            mw_difference = None
            mw_valid = True

            if target_mw:
                try:
                    target_mw_float = float(target_mw)
                    mw_difference = round(abs(calculated_mw - target_mw_float), 2)
                    mw_valid = mw_difference <= 1.0
                except (ValueError, TypeError):
                    pass

            return {
                "status": "valid" if mw_valid else "mw_mismatch",
                "desalted_smiles": desalted_smiles,
                "calculated_mw": calculated_mw,
                "mw_difference": mw_difference,
                "mw_valid": mw_valid
            }

        except ImportError:
            # RDKit not installed - return basic validation
            logger.warning("RDKit not installed. Skipping MW validation.")
            return {
                "status": "unverified",
                "desalted_smiles": smiles,
                "calculated_mw": None,
                "mw_difference": None,
                "note": "RDKit not available for validation"
            }
        except Exception as e:
            logger.error(f"RDKit validation error: {e}")
            return {
                "status": "error",
                "error": str(e)
            }

    async def autofill_smiles(self, record: Dict[str, Any]) -> Dict[str, Any]:
        """
        레코드의 CAS 또는 Name으로 SMILES 자동 채우기
        우선순위: CAS Number > Product Name
        """
        cas_number = record.get("cas_number")
        product_name = record.get("product_name") or record.get("name")
        target_mw = record.get("molecular_weight")

        result = None

        # Try CAS first
        if cas_number:
            result = await self.get_smiles_by_cas(cas_number)
            if "error" not in result:
                # Additional MW validation against record's MW
                if target_mw and result.get("calculated_mw"):
                    try:
                        diff = abs(float(result["calculated_mw"]) - float(target_mw))
                        result["target_mw_match"] = diff <= 1.0
                        result["target_mw_diff"] = round(diff, 2)
                    except:
                        pass
                return result

        # Fallback to name
        if product_name:
            result = await self.get_smiles_by_name(product_name)
            if "error" not in result:
                if target_mw and result.get("calculated_mw"):
                    try:
                        diff = abs(float(result["calculated_mw"]) - float(target_mw))
                        result["target_mw_match"] = diff <= 1.0
                        result["target_mw_diff"] = round(diff, 2)
                    except:
                        pass
                return result

        return {"error": "No CAS number or product name available for lookup"}

    def clear_cache(self):
        """캐시 클리어"""
        self._cache.clear()
        logger.info("PubChem cache cleared")


# Singleton instance
pubchem_service = PubChemService()
