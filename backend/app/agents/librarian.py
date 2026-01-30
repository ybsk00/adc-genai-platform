"""
The Librarian Agent - Dynamic Substructure RAG
스캐폴드 추출 및 지식베이스 의미론적 매핑

Phase 1 Enhancement:
- Gemini 3.0 Pro 통합 (evidence synthesis)
- Real-time Streaming UI 연동
- Target Normalization 활용
"""
from typing import List, Dict, Any, Optional
import logging
import json
from datetime import datetime

from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_core.messages import HumanMessage, SystemMessage

from app.agents.base_agent import BaseDesignAgent, AgentOutput
from app.agents.design_state import DesignSessionState
from app.services.rag_service import RAGService
from app.core.supabase import get_supabase_client
from app.core.config import settings
from app.core.websocket_hub import websocket_hub

logger = logging.getLogger(__name__)


class LibrarianAgent(BaseDesignAgent):
    """
    The Librarian: Dynamic Substructure RAG 에이전트

    핵심 기능:
    1. 입력 SMILES에서 스캐폴드/서브스트럭처 추출 (Murcko decomposition)
    2. 스캐폴드를 벡터 임베딩하여 의미론적 유사성 검색
    3. knowledge_base + golden_set에서 관련 문헌/구조 검색
    4. 관련 PMID 및 근거 자료 제공
    5. [Phase 1] Gemini를 활용한 근거 합성 및 설명 생성

    Phase 1 Enhancement:
    - Gemini Pro를 사용한 evidence synthesis
    - Target Normalization을 통한 정확한 검색
    - Real-time Streaming을 통한 진행 상황 전송
    """

    name = "librarian"

    def __init__(self):
        super().__init__()
        self.rag_service = RAGService()
        # Gemini 모델 초기화 (Phase 1)
        self.llm = ChatGoogleGenerativeAI(
            model=settings.GEMINI_PRO_MODEL_ID,  # gemini-2.0-flash
            google_api_key=settings.GOOGLE_API_KEY,
            temperature=0.2  # 낮은 temperature로 일관된 근거 생성
        )

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
            # [Streaming] 시작 로그
            await websocket_hub.stream_agent_log(
                session_id, "info",
                f"문헌 검색 시작: 타겟 '{target_antigen}'",
                emoji="📚", agent_name="librarian"
            )

            # 1. 스캐폴드 추출
            await websocket_hub.stream_progress(session_id, 10, "스캐폴드 추출 중")
            scaffolds = self._extract_scaffolds(smiles)
            await websocket_hub.stream_agent_log(
                session_id, "success",
                f"{len(scaffolds)}개 스캐폴드 추출 완료",
                emoji="✅", agent_name="librarian"
            )

            # 2. 타겟 정규화 활용 (Phase 1)
            normalized_target = await self._normalize_target(target_antigen)
            if normalized_target != target_antigen:
                await websocket_hub.stream_agent_log(
                    session_id, "info",
                    f"타겟 정규화: '{target_antigen}' → '{normalized_target}'",
                    emoji="🔄", agent_name="librarian"
                )

            # 3. 각 스캐폴드에 대해 지식베이스 매핑
            await websocket_hub.stream_progress(session_id, 30, "지식베이스 검색 중")
            knowledge_results = []
            for idx, scaffold in enumerate(scaffolds):
                # 3.1 기존 매핑 검색 또는 새로 생성
                mapping = await self._get_or_create_scaffold_mapping(scaffold)

                # 3.2 의미론적 검색
                semantic_hits = await self._semantic_search(
                    scaffold["smiles"],
                    normalized_target,
                    top_k=5
                )

                # 3.3 Golden Set 유사도 기반 검색
                golden_set_hits = await self._search_golden_set_by_similarity(
                    normalized_target,
                    threshold=0.6
                )

                knowledge_results.append({
                    "scaffold_type": scaffold["type"],
                    "scaffold_smiles": scaffold["smiles"],
                    "semantic_hits": semantic_hits,
                    "golden_set_hits": golden_set_hits,
                    "mapping_id": str(mapping["id"]) if mapping else None
                })

                await websocket_hub.stream_progress(
                    session_id,
                    30 + (idx + 1) / len(scaffolds) * 30,
                    f"스캐폴드 {idx + 1}/{len(scaffolds)} 검색 완료"
                )

            # 4. 통합 근거 생성
            await websocket_hub.stream_progress(session_id, 70, "근거 자료 통합 중")
            evidence = self._compile_evidence(knowledge_results)

            # 5. [Phase 1] Gemini를 사용한 근거 합성
            await websocket_hub.stream_progress(session_id, 80, "AI 근거 합성 중")
            await websocket_hub.stream_agent_log(
                session_id, "info",
                "Gemini Pro로 근거 합성 중...",
                emoji="🤖", agent_name="librarian"
            )
            synthesized_evidence = await self._synthesize_evidence_with_gemini(
                target_antigen=normalized_target,
                scaffolds=scaffolds,
                semantic_hits=[h for kr in knowledge_results for h in kr.get("semantic_hits", [])],
                golden_set_hits=[h for kr in knowledge_results for h in kr.get("golden_set_hits", [])]
            )

            await websocket_hub.stream_progress(session_id, 100, "검색 완료")
            await websocket_hub.stream_agent_log(
                session_id, "success",
                f"검색 완료: {len(evidence['pmids'])} 문헌, {len(evidence['golden_sets'])} Golden Set 참조",
                emoji="✅", agent_name="librarian"
            )

            # [Streaming] 참조 문헌 전송
            await websocket_hub.broadcast_librarian_references(
                session_id,
                references=[{
                    "id": h.get("id"),
                    "pmid": h.get("pmid"),
                    "title": h.get("title", ""),
                    "relevance_score": h.get("relevance_score", 0),
                    "summary": h.get("summary", "")
                } for kr in knowledge_results for h in kr.get("semantic_hits", [])],
                golden_set_refs=evidence["golden_sets"]
            )

            output = AgentOutput(
                success=True,
                data={
                    "scaffolds": scaffolds,
                    "knowledge_results": knowledge_results,
                    "evidence_summary": synthesized_evidence.get("summary", evidence["summary"]),
                    "evidence_reasoning": synthesized_evidence.get("reasoning", ""),
                    "pmid_references": evidence["pmids"],
                    "golden_set_references": evidence["golden_sets"],
                    "normalized_target": normalized_target
                },
                reasoning=synthesized_evidence.get("reasoning",
                    f"Extracted {len(scaffolds)} scaffolds, found {len(evidence['pmids'])} literature references"),
                confidence_score=0.85 if evidence["pmids"] else 0.6,
                referenced_knowledge_ids=evidence["knowledge_ids"],
                referenced_pmids=evidence["pmids"]
            )

            await self._log_complete(session_id, output)
            return output

        except Exception as e:
            logger.exception(f"[librarian] Error: {e}")
            await websocket_hub.stream_agent_log(
                session_id, "error",
                f"검색 오류: {str(e)}",
                emoji="❌", agent_name="librarian"
            )
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

    # ============== Phase 1: Gemini Integration ==============

    async def _normalize_target(self, target: str) -> str:
        """
        [Phase 1] Target Normalization
        target_synonyms 테이블을 사용하여 타겟명을 정규화합니다.
        """
        if not target:
            return target

        try:
            # RPC 함수 호출 (normalize_target)
            result = self.supabase.rpc("normalize_target", {"input_target": target}).execute()
            if result.data:
                return result.data
            return target
        except Exception as e:
            logger.warning(f"[librarian] Target normalization failed: {e}")
            return target

    async def _synthesize_evidence_with_gemini(
        self,
        target_antigen: str,
        scaffolds: List[Dict],
        semantic_hits: List[Dict],
        golden_set_hits: List[Dict]
    ) -> Dict[str, Any]:
        """
        [Phase 1] Gemini를 사용하여 수집된 근거를 합성합니다.

        검색 결과를 바탕으로 구조화된 근거 요약을 생성합니다.
        """
        if not (semantic_hits or golden_set_hits):
            return {
                "summary": "관련 문헌을 찾지 못했습니다.",
                "reasoning": "검색 결과가 없어 근거 합성을 수행할 수 없습니다."
            }

        system_prompt = """당신은 ADC(Antibody-Drug Conjugate) 개발 전문가입니다.
검색된 문헌과 Golden Set 참조를 분석하여 설계 근거를 합성해야 합니다.

응답 형식 (JSON):
{
    "summary": "핵심 근거 요약 (2-3문장)",
    "reasoning": "상세 분석 내용",
    "key_findings": ["주요 발견사항 리스트"],
    "recommendations": ["설계 권장사항 리스트"],
    "confidence": 0.0~1.0 (근거 신뢰도)
}"""

        # 검색 결과 요약
        literature_summary = "\n".join([
            f"- [{h.get('pmid', 'N/A')}] {h.get('title', 'Unknown')}: {h.get('summary', '')[:100]}..."
            for h in semantic_hits[:5]
        ]) if semantic_hits else "검색된 문헌 없음"

        golden_set_summary = "\n".join([
            f"- {h.get('name', 'Unknown')}: Target={h.get('target', 'N/A')}, Status={h.get('status', 'N/A')}"
            for h in golden_set_hits[:5]
        ]) if golden_set_hits else "검색된 Golden Set 없음"

        scaffold_summary = "\n".join([
            f"- {s.get('type', 'unknown')}: {s.get('smiles', '')[:50]}..."
            for s in scaffolds[:3]
        ]) if scaffolds else "추출된 스캐폴드 없음"

        user_prompt = f"""타겟 항원: {target_antigen}

추출된 스캐폴드:
{scaffold_summary}

관련 문헌:
{literature_summary}

Golden Set 참조 (FDA 승인 ADC):
{golden_set_summary}

위 정보를 바탕으로 ADC 설계 근거를 합성해주세요."""

        try:
            response = await self.llm.ainvoke([
                SystemMessage(content=system_prompt),
                HumanMessage(content=user_prompt)
            ])

            content = response.content
            # JSON 파싱 시도
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0]
            elif "```" in content:
                content = content.split("```")[1].split("```")[0]

            data = json.loads(content)
            return {
                "summary": data.get("summary", ""),
                "reasoning": data.get("reasoning", ""),
                "key_findings": data.get("key_findings", []),
                "recommendations": data.get("recommendations", []),
                "confidence": data.get("confidence", 0.7)
            }

        except json.JSONDecodeError as e:
            logger.warning(f"[librarian] Gemini JSON parse error: {e}")
            # Fallback: 원본 응답 사용
            return {
                "summary": response.content[:500] if response else "",
                "reasoning": "Gemini 응답 파싱 실패"
            }

        except Exception as e:
            logger.error(f"[librarian] Gemini synthesis error: {e}")
            return {
                "summary": f"문헌 {len(semantic_hits)}개, Golden Set {len(golden_set_hits)}개 검색됨",
                "reasoning": f"AI 합성 중 오류 발생: {str(e)}"
            }

    async def _search_antibody_library(
        self,
        target: str,
        top_k: int = 5
    ) -> List[Dict]:
        """
        [Phase 1] antibody_library에서 타겟 기반 항체 검색

        target_normalized 컬럼을 활용한 검색
        """
        try:
            result = self.supabase.table("antibody_library").select(
                "id, product_name, cat_no, related_disease, target_normalized"
            ).or_(
                f"target_normalized.ilike.%{target}%,related_disease.ilike.%{target}%"
            ).limit(top_k).execute()

            return result.data or []

        except Exception as e:
            logger.warning(f"[librarian] Antibody library search error: {e}")
            return []

    # ============== One-Click Navigator Support ==============

    async def find_antibodies_by_disease(
        self,
        disease_name: str,
        top_k: int = 3
    ) -> List[Dict[str, Any]]:
        """
        [One-Click Navigator] 질환명 기반 최적 항체 검색

        Args:
            disease_name: 질환명 (예: "Pancreatic Cancer")
            top_k: 반환할 항체 수

        Returns:
            List[Dict]: 항체 후보 목록
        """
        try:
            # 1. 질환명 임베딩 생성
            disease_embedding = await self.rag_service.generate_embedding(
                f"Disease: {disease_name}, treatment target proteins, therapeutic antibodies"
            )

            # 2. antibody_library 벡터 검색
            results = self.supabase.rpc("match_antibody_by_disease", {
                "query_embedding": disease_embedding,
                "match_threshold": 0.5,
                "match_count": top_k * 2
            }).execute()

            antibodies = results.data or []

            # 3. 임상 데이터 기반 재순위화
            scored_antibodies = []
            for ab in antibodies:
                clinical_score = self._calculate_antibody_clinical_score(ab)
                scored_antibodies.append({
                    **ab,
                    "clinical_score": clinical_score,
                    "combined_score": ab.get("similarity", 0.5) * 0.4 + clinical_score * 0.6
                })

            # 4. Top K 선정
            top_antibodies = sorted(
                scored_antibodies,
                key=lambda x: x["combined_score"],
                reverse=True
            )[:top_k]

            # 5. Fallback: 벡터 검색 결과 없으면 직접 검색
            if not top_antibodies:
                logger.info(f"[librarian] Vector search empty, trying direct search: {disease_name}")
                return await self._direct_disease_antibody_search(disease_name, top_k)

            return top_antibodies

        except Exception as e:
            logger.error(f"[librarian] find_antibodies_by_disease error: {e}")
            return await self._direct_disease_antibody_search(disease_name, top_k)

    async def _direct_disease_antibody_search(
        self,
        disease_name: str,
        top_k: int
    ) -> List[Dict[str, Any]]:
        """직접 키워드 검색 (fallback)"""
        try:
            result = self.supabase.table("antibody_library").select(
                "id, product_name, target_normalized, isotype, related_disease, full_spec"
            ).ilike(
                "related_disease", f"%{disease_name}%"
            ).limit(top_k).execute()

            if result.data:
                return [
                    {
                        "id": ab["id"],
                        "name": ab["product_name"],
                        "target_protein": ab.get("target_normalized", "Unknown"),
                        "isotype": ab.get("isotype"),
                        "related_disease": ab.get("related_disease"),
                        "full_spec": ab.get("full_spec"),
                        "clinical_score": 0.5,
                        "combined_score": 0.5,
                        "similarity": 0.5
                    }
                    for ab in result.data
                ]

            # Golden Set에서 검색
            gs_result = self.supabase.table("golden_set").select(
                "id, drug_name, target_1, indication, orr_pct, os_months"
            ).ilike(
                "indication", f"%{disease_name}%"
            ).limit(top_k).execute()

            return [
                {
                    "id": gs["id"],
                    "name": gs["drug_name"],
                    "target_protein": gs.get("target_1", "Unknown"),
                    "isotype": None,
                    "related_disease": gs.get("indication"),
                    "full_spec": None,
                    "orr_pct": gs.get("orr_pct"),
                    "os_months": gs.get("os_months"),
                    "clinical_score": 0.7,
                    "combined_score": 0.7,
                    "similarity": 0.7
                }
                for gs in (gs_result.data or [])
            ]

        except Exception as e:
            logger.error(f"[librarian] Direct disease search error: {e}")
            return []

    def _calculate_antibody_clinical_score(self, ab: Dict[str, Any]) -> float:
        """항체의 임상 점수 계산"""
        score = 0.0

        # ORR (40%)
        orr = ab.get("orr_pct", 0) or 0
        score += (orr / 100) * 0.4

        # OS (30%)
        os_months = ab.get("os_months", 0) or 0
        score += min(os_months / 36, 1.0) * 0.3

        # 기본 점수 (30%)
        score += 0.3

        return min(score, 1.0)
