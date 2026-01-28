"""
The Librarian Agent - Dynamic Substructure RAG
스캐폴드 추출 및 지식베이스 의미론적 매핑
"""
from typing import List, Dict, Any, Optional
import logging
from datetime import datetime

from app.agents.base_agent import BaseDesignAgent, AgentOutput
from app.agents.design_state import DesignSessionState
from app.services.rag_service import RAGService
from app.core.supabase import get_supabase_client

logger = logging.getLogger(__name__)


class LibrarianAgent(BaseDesignAgent):
    """
    The Librarian: Dynamic Substructure RAG 에이전트

    핵심 기능:
    1. 입력 SMILES에서 스캐폴드/서브스트럭처 추출 (Murcko decomposition)
    2. 스캐폴드를 벡터 임베딩하여 의미론적 유사성 검색
    3. knowledge_base + golden_set에서 관련 문헌/구조 검색
    4. 관련 PMID 및 근거 자료 제공
    """

    name = "librarian"

    def __init__(self):
        super().__init__()
        self.rag_service = RAGService()

    async def execute(self, state: DesignSessionState) -> AgentOutput:
        """메인 실행"""
        session_id = state["session_id"]
        smiles = state["current_smiles"]
        target_antigen = state.get("target_antigen", "")

        await self._log_start(session_id, {
            "smiles_length": len(smiles) if smiles else 0,
            "target": target_antigen
        })

        try:
            # 1. 스캐폴드 추출
            scaffolds = self._extract_scaffolds(smiles)

            # 2. 각 스캐폴드에 대해 지식베이스 매핑
            knowledge_results = []
            for scaffold in scaffolds:
                # 2.1 기존 매핑 검색 또는 새로 생성
                mapping = await self._get_or_create_scaffold_mapping(scaffold)

                # 2.2 의미론적 검색
                semantic_hits = await self._semantic_search(
                    scaffold["smiles"],
                    target_antigen,
                    top_k=5
                )

                # 2.3 Golden Set 유사도 기반 검색
                golden_set_hits = await self._search_golden_set_by_similarity(
                    target_antigen,
                    threshold=0.6
                )

                knowledge_results.append({
                    "scaffold_type": scaffold["type"],
                    "scaffold_smiles": scaffold["smiles"],
                    "semantic_hits": semantic_hits,
                    "golden_set_hits": golden_set_hits,
                    "mapping_id": str(mapping["id"]) if mapping else None
                })

            # 3. 통합 근거 생성
            evidence = self._compile_evidence(knowledge_results)

            output = AgentOutput(
                success=True,
                data={
                    "scaffolds": scaffolds,
                    "knowledge_results": knowledge_results,
                    "evidence_summary": evidence["summary"],
                    "pmid_references": evidence["pmids"],
                    "golden_set_references": evidence["golden_sets"]
                },
                reasoning=f"Extracted {len(scaffolds)} scaffolds, found {len(evidence['pmids'])} literature references",
                confidence_score=0.8 if evidence["pmids"] else 0.5,
                referenced_knowledge_ids=evidence["knowledge_ids"],
                referenced_pmids=evidence["pmids"]
            )

            await self._log_complete(session_id, output)
            return output

        except Exception as e:
            logger.exception(f"[librarian] Error: {e}")
            await self._log_error(session_id, str(e))
            return AgentOutput(
                success=False,
                data={},
                error=str(e)
            )

    def _extract_scaffolds(self, smiles: str) -> List[Dict[str, Any]]:
        """
        스캐폴드 추출 (RDKit Murcko Decomposition)

        Returns: [{type, smiles}, ...]
        """
        if not smiles:
            return []

        scaffolds = []

        try:
            from rdkit import Chem
            from rdkit.Chem.Scaffolds import MurckoScaffold

            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                logger.warning(f"[librarian] Invalid SMILES: {smiles[:50]}")
                return []

            # 1. Generic Scaffold
            try:
                generic = MurckoScaffold.MakeScaffoldGeneric(
                    MurckoScaffold.GetScaffoldForMol(mol)
                )
                if generic:
                    generic_smiles = Chem.MolToSmiles(generic)
                    scaffolds.append({
                        "type": "murcko_generic",
                        "smiles": generic_smiles
                    })
            except Exception as e:
                logger.warning(f"[librarian] Generic scaffold error: {e}")

            # 2. Framework
            try:
                framework = MurckoScaffold.GetScaffoldForMol(mol)
                if framework:
                    framework_smiles = Chem.MolToSmiles(framework)
                    scaffolds.append({
                        "type": "murcko_framework",
                        "smiles": framework_smiles
                    })
            except Exception as e:
                logger.warning(f"[librarian] Framework error: {e}")

            # 3. Ring Systems (간단한 추출)
            try:
                ring_info = mol.GetRingInfo()
                if ring_info.NumRings() > 0:
                    scaffolds.append({
                        "type": "ring_system",
                        "smiles": smiles,  # 원본 유지
                        "ring_count": ring_info.NumRings()
                    })
            except Exception as e:
                logger.warning(f"[librarian] Ring system error: {e}")

        except ImportError:
            logger.warning("[librarian] RDKit not available, skipping scaffold extraction")
            # Fallback: 원본 SMILES 반환
            scaffolds.append({
                "type": "original",
                "smiles": smiles
            })

        except Exception as e:
            logger.error(f"[librarian] Scaffold extraction error: {e}")

        return scaffolds

    async def _get_or_create_scaffold_mapping(self, scaffold: Dict) -> Optional[Dict]:
        """스캐폴드 매핑 조회 또는 신규 생성"""
        try:
            # 기존 매핑 검색
            result = self.supabase.table("scaffold_knowledge_mapping").select("*").eq(
                "scaffold_smiles", scaffold["smiles"]
            ).eq(
                "scaffold_type", scaffold["type"]
            ).single().execute()

            if result.data:
                # 조회 카운트 증가
                self.supabase.table("scaffold_knowledge_mapping").update({
                    "retrieval_count": result.data["retrieval_count"] + 1,
                    "last_retrieved_at": datetime.utcnow().isoformat()
                }).eq("id", result.data["id"]).execute()
                return result.data

            # 신규 매핑 생성
            new_mapping = self.supabase.table("scaffold_knowledge_mapping").insert({
                "scaffold_smiles": scaffold["smiles"],
                "scaffold_type": scaffold["type"],
                "scaffold_description": f"Scaffold of type {scaffold['type']}",
                "retrieval_count": 1,
                "last_retrieved_at": datetime.utcnow().isoformat()
            }).execute()

            return new_mapping.data[0] if new_mapping.data else None

        except Exception as e:
            logger.warning(f"[librarian] Scaffold mapping error: {e}")
            return None

    async def _semantic_search(
        self,
        scaffold_smiles: str,
        target: str,
        top_k: int = 5
    ) -> List[Dict]:
        """의미론적 검색: 스캐폴드 + 타겟 기반 knowledge_base 검색"""
        query = f"ADC scaffold structure {scaffold_smiles[:50]} targeting {target}"

        try:
            # RAG 검색
            results = await self.rag_service.search(
                query=query,
                top_k=top_k
            )

            return [{
                "id": r.get("id"),
                "title": r.get("title", ""),
                "relevance_score": r.get("score", 0),
                "pmid": r.get("pmid") or r.get("metadata", {}).get("pmid"),
                "summary": r.get("summary", "")[:200]
            } for r in results]

        except Exception as e:
            logger.warning(f"[librarian] Semantic search error: {e}")
            return []

    async def _search_golden_set_by_similarity(
        self,
        target: str,
        threshold: float = 0.6
    ) -> List[Dict]:
        """Golden Set에서 타겟 기반 검색"""
        try:
            result = self.supabase.table("golden_set").select(
                "id, drug_name, target_1, target_2, indication, clinical_status"
            ).or_(
                f"target_1.ilike.%{target}%,target_2.ilike.%{target}%"
            ).limit(10).execute()

            return [{
                "id": str(g["id"]),
                "name": g.get("drug_name"),
                "target": g.get("target_1"),
                "status": g.get("clinical_status")
            } for g in (result.data or [])]

        except Exception as e:
            logger.warning(f"[librarian] Golden Set search error: {e}")
            return []

    def _compile_evidence(self, knowledge_results: List[Dict]) -> Dict[str, Any]:
        """검색 결과를 통합하여 근거 자료 생성"""
        pmids = set()
        knowledge_ids = []
        golden_sets = []
        summaries = []

        for result in knowledge_results:
            for hit in result.get("semantic_hits", []):
                if hit.get("pmid"):
                    pmids.add(hit["pmid"])
                if hit.get("id"):
                    knowledge_ids.append(hit["id"])

            for gs in result.get("golden_set_hits", []):
                golden_sets.append({
                    "id": gs["id"],
                    "name": gs.get("name"),
                    "target": gs.get("target")
                })

            summaries.append(
                f"{result['scaffold_type']}: {len(result.get('semantic_hits', []))} refs"
            )

        return {
            "pmids": list(pmids),
            "knowledge_ids": knowledge_ids,
            "golden_sets": golden_sets[:5],
            "summary": "; ".join(summaries) if summaries else "No references found"
        }

    async def link_scaffold_to_knowledge(
        self,
        scaffold_mapping_id: str,
        knowledge_ids: List[int],
        golden_set_ids: List[str],
        pmids: List[str]
    ):
        """스캐폴드 매핑에 지식베이스 링크 추가"""
        try:
            self.supabase.table("scaffold_knowledge_mapping").update({
                "linked_knowledge_ids": knowledge_ids,
                "linked_golden_set_ids": golden_set_ids,
                "linked_pmids": pmids,
                "updated_at": datetime.utcnow().isoformat()
            }).eq("id", scaffold_mapping_id).execute()
        except Exception as e:
            logger.error(f"[librarian] Failed to link scaffold: {e}")
