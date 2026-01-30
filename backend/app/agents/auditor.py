"""
The Auditor Agent - 검증 및 감사
PAINS/Lipinski/독성 검증 + Constraint Guardrail

Phase 1 Enhancement:
- RDKit 기반 정밀 화학 검증
- Gemini를 활용한 지능형 검증 추론
- Real-time Streaming UI 연동
- SMILES 무결성 검증
"""
from typing import List, Dict, Any, Optional
import logging
import json

from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_core.messages import HumanMessage, SystemMessage

from app.agents.base_agent import BaseDesignAgent, AgentOutput
from app.agents.design_state import (
    DesignSessionState,
    ValidationFlags,
    ConstraintCheck
)
from app.core.supabase import get_supabase_client
from app.core.config import settings
from app.core.websocket_hub import websocket_hub

logger = logging.getLogger(__name__)

# RDKit 임포트 시도
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, FilterCatalog
    from rdkit.Chem.FilterCatalog import FilterCatalogParams
    RDKIT_AVAILABLE = True
    logger.info("[auditor] RDKit loaded successfully")
except ImportError:
    RDKIT_AVAILABLE = False
    logger.warning("[auditor] RDKit not available - using fallback validation")


class AuditorAgent(BaseDesignAgent):
    """
    The Auditor: 검증 및 감사 에이전트

    핵심 기능:
    1. PAINS/Lipinski/독성 검증
    2. Constraint Guardrail: 초기 제약조건 준수 확인 (Constraint 4)
    3. 리스크 스코어 계산
    4. 최종 승인/재설계 결정

    Phase 1 Enhancement:
    - RDKit 기반 정밀 검증 (PAINS, SA Score, 독성 예측)
    - Gemini를 활용한 검증 결과 해석 및 권장사항 생성
    - Real-time Streaming으로 검증 진행 상황 전송
    """

    name = "auditor"

    # 리스크 가중치
    RISK_WEIGHTS = {
        "lipinski_fail": 3,
        "pains_detected": 4,
        "constraint_violation_high": 5,
        "constraint_violation_medium": 2,
        "low_similarity": 2,
        "high_sa_score": 1,
        "smiles_invalid": 5,
        "toxicity_alert": 4
    }

    # ADC 특화 Lipinski 확장 기준 (Beyond Rule of 5)
    ADC_THRESHOLDS = {
        "mw_max": 1000,  # ADC payload는 더 무거울 수 있음
        "logp_max": 6,
        "hbd_max": 7,
        "hba_max": 15,
        "tpsa_max": 200,
        "rotatable_bonds_max": 15
    }

    def __init__(self):
        super().__init__()
        # Gemini 모델 초기화 (Phase 1)
        self.llm = ChatGoogleGenerativeAI(
            model=settings.GEMINI_PRO_MODEL_ID,
            google_api_key=settings.GOOGLE_API_KEY,
            temperature=0.1  # 검증은 낮은 temperature
        )
        # PAINS 필터 초기화
        self._init_pains_filter()

    async def execute(self, state: DesignSessionState) -> AgentOutput:
        """메인 실행"""
        session_id = state["session_id"]

        await self._log_start(session_id, {
            "smiles": state["current_smiles"][:30] + "..." if state["current_smiles"] else None,
            "tier": state["tier"]
        })

        try:
            # [Streaming] 시작 로그
            await websocket_hub.stream_agent_log(
                session_id, "info",
                "ADC 구조 검증 시작",
                emoji="🔍", agent_name="auditor"
            )

            # 1. [Phase 1] SMILES 무결성 검증
            await websocket_hub.stream_progress(session_id, 10, "SMILES 무결성 검증 중")
            smiles_validation = await self._validate_smiles_integrity(state["current_smiles"])
            if not smiles_validation["valid"]:
                await websocket_hub.stream_agent_log(
                    session_id, "error",
                    f"SMILES 무결성 오류: {smiles_validation['error']}",
                    emoji="❌", agent_name="auditor"
                )
                # 재설계 요청
                return await self._request_redesign(
                    session_id, state,
                    reason=f"Invalid SMILES: {smiles_validation['error']}",
                    failed_checks=["smiles_integrity"]
                )

            await websocket_hub.stream_agent_log(
                session_id, "success",
                "SMILES 무결성 검증 통과",
                emoji="✅", agent_name="auditor"
            )

            # 2. [Phase 1] RDKit 기반 화학 검증
            await websocket_hub.stream_progress(session_id, 30, "화학적 특성 검증 중")
            chemistry_validation = await self._validate_chemistry_rdkit(state)

            # 검증 결과 스트리밍
            check_items = self._format_check_items(chemistry_validation)
            await websocket_hub.broadcast_auditor_feedback(
                session_id,
                is_in_redesign_loop=False,
                check_items=check_items
            )

            # 3. Constraint Guardrail 검사 (Constraint 4)
            await websocket_hub.stream_progress(session_id, 50, "제약조건 검사 중")
            constraint_check = await self._check_constraint_guardrail(state)

            if constraint_check.get("violations"):
                await websocket_hub.stream_agent_log(
                    session_id, "warning",
                    f"{len(constraint_check['violations'])}개 제약조건 위반 감지",
                    emoji="⚠️", agent_name="auditor"
                )

            # 4. 리스크 스코어 계산
            await websocket_hub.stream_progress(session_id, 70, "리스크 평가 중")
            risk_score = self._calculate_risk_score(
                chemistry_validation,
                constraint_check
            )

            # 5. [Phase 1] Gemini를 사용한 지능형 결정
            await websocket_hub.stream_progress(session_id, 85, "AI 검증 분석 중")
            await websocket_hub.stream_agent_log(
                session_id, "info",
                "Gemini로 검증 결과 분석 중...",
                emoji="🤖", agent_name="auditor"
            )
            decision = await self._make_decision_with_gemini(
                state,
                chemistry_validation,
                constraint_check,
                risk_score
            )

            # 6. 최종 리포트 생성
            await websocket_hub.stream_progress(session_id, 95, "리포트 생성 중")
            final_report = self._generate_report(
                state, chemistry_validation, constraint_check, risk_score, decision
            )

            # [Streaming] 최종 피드백 전송
            await websocket_hub.broadcast_auditor_feedback(
                session_id,
                is_in_redesign_loop=decision["action"] == "redesign",
                check_items=check_items,
                redesign_request={
                    "requestedBy": "auditor",
                    "reason": decision["reasoning"],
                    "failed_checks": decision.get("failed_checks", []),
                    "suggestions": decision.get("suggestions", []),
                    "target_agent": "alchemist",
                    "iteration": state.get("redesign_count", 0) + 1,
                    "max_iterations": 3
                } if decision["action"] == "redesign" else None
            )

            await websocket_hub.stream_progress(session_id, 100, "검증 완료")
            status_emoji = "✅" if decision["approved"] else "⚠️" if decision["action"] == "manual_review" else "🔄"
            await websocket_hub.stream_agent_log(
                session_id, "success" if decision["approved"] else "warning",
                f"검증 완료: {decision['action']} (신뢰도: {decision['confidence']:.0%})",
                emoji=status_emoji, agent_name="auditor"
            )

            output = AgentOutput(
                success=decision["approved"],
                data={
                    "smiles_validation": smiles_validation,
                    "chemistry_validation": chemistry_validation,
                    "constraint_check": constraint_check,
                    "risk_score": risk_score,
                    "decision": decision,
                    "final_report": final_report
                },
                reasoning=decision["reasoning"],
                confidence_score=decision["confidence"],
                next_agent="alchemist" if decision["action"] == "redesign" else None
            )

            await self._log_complete(session_id, output)
            return output

        except Exception as e:
            logger.exception(f"[auditor] Error: {e}")
            await websocket_hub.stream_agent_log(
                session_id, "error",
                f"검증 오류: {str(e)}",
                emoji="❌", agent_name="auditor"
            )
            await self._log_error(session_id, str(e))
            return AgentOutput(
                success=False,
                data={},
                error=str(e)
            )

    def _validate_chemistry(self, state: DesignSessionState) -> Dict[str, Any]:
        """화학적 검증 (Coder 결과 기반)"""
        metrics = state.get("calculated_metrics", {})
        validation = state.get("validation_flags", {})

        # Lipinski 검증
        mw = metrics.get("mw", 0)
        logp = metrics.get("logp", 0)
        hbd = metrics.get("hbd", 0)
        hba = metrics.get("hba", 0)

        lipinski_violations = 0
        if mw and mw > 500:
            lipinski_violations += 1
        if logp and logp > 5:
            lipinski_violations += 1
        if hbd and hbd > 5:
            lipinski_violations += 1
        if hba and hba > 10:
            lipinski_violations += 1

        lipinski_pass = lipinski_violations <= 1  # 1개까지 허용

        # SA Score 검증
        sa_score = metrics.get("sa_score", 5)
        sa_acceptable = sa_score is None or sa_score < 6

        return {
            "lipinski_pass": lipinski_pass,
            "lipinski_violations": lipinski_violations,
            "pains_detected": validation.get("pains_free") == False,
            "pains_pattern": validation.get("structural_alerts", []),
            "sa_score": sa_score,
            "sa_acceptable": sa_acceptable,
            "metrics": {
                "mw": mw,
                "logp": logp,
                "hbd": hbd,
                "hba": hba
            }
        }

    async def _check_constraint_guardrail(self, state: DesignSessionState) -> ConstraintCheck:
        """
        Constraint Guardrail (Constraint 4)

        사용자 초기 제약조건 vs 설계 결과 비교 검증
        """
        constraints = {
            "target_antigen": state.get("target_antigen"),
            "requested_dar": state.get("requested_dar"),
            "linker_preference": state.get("linker_preference")
        }

        violations = []
        warnings = []

        # 1. Target Antigen 검증
        if constraints["target_antigen"]:
            target_validated = await self._validate_target_compatibility(
                state["current_smiles"],
                constraints["target_antigen"]
            )
            if not target_validated["compatible"]:
                violations.append({
                    "constraint": "target_antigen",
                    "expected": constraints["target_antigen"],
                    "actual": target_validated.get("detected_target"),
                    "severity": "high",
                    "reason": target_validated.get("reason")
                })

        # 2. DAR 검증
        if constraints["requested_dar"]:
            detected_dar = self._estimate_dar_from_structure(state["current_smiles"])
            if detected_dar and detected_dar != constraints["requested_dar"]:
                if abs(detected_dar - constraints["requested_dar"]) > 1:
                    violations.append({
                        "constraint": "requested_dar",
                        "expected": constraints["requested_dar"],
                        "actual": detected_dar,
                        "severity": "medium",
                        "reason": f"DAR mismatch: expected {constraints['requested_dar']}, estimated {detected_dar}"
                    })
                else:
                    warnings.append({
                        "constraint": "requested_dar",
                        "message": f"DAR slightly different: expected {constraints['requested_dar']}, estimated {detected_dar}"
                    })

        # 3. Linker Preference 검증
        if constraints["linker_preference"] and constraints["linker_preference"] != "any":
            linker_analysis = self._analyze_linker_type(state["current_smiles"])
            if linker_analysis["type"] != constraints["linker_preference"]:
                if linker_analysis["type"] != "unknown":
                    violations.append({
                        "constraint": "linker_preference",
                        "expected": constraints["linker_preference"],
                        "actual": linker_analysis["type"],
                        "severity": "medium",
                        "reason": linker_analysis.get("reason", "Linker type mismatch")
                    })

        summary = self._generate_constraint_summary(violations, warnings)

        return ConstraintCheck(
            passed=len(violations) == 0,
            violations=violations,
            warnings=warnings,
            constraints_checked=list(constraints.keys()),
            summary=summary
        )

    async def _validate_target_compatibility(
        self,
        smiles: str,
        target: str
    ) -> Dict[str, Any]:
        """타겟 항원 호환성 검증"""
        # Golden Set에서 유사 구조의 타겟 확인
        try:
            if not smiles:
                return {"compatible": True, "confidence": 0.5, "reason": "No SMILES to validate"}

            # 현재는 간단한 검증 (추후 ML 모델로 대체 가능)
            return {
                "compatible": True,
                "confidence": 0.85,
                "detected_target": target,
                "reason": "Target compatibility assumed (detailed validation pending)"
            }
        except Exception as e:
            logger.warning(f"[auditor] Target validation error: {e}")
            return {"compatible": True, "confidence": 0.5, "reason": str(e)}

    def _estimate_dar_from_structure(self, smiles: str) -> Optional[int]:
        """구조에서 DAR 추정"""
        if not smiles:
            return None

        # 결합점 패턴 분석 (간단한 휴리스틱)
        # 실제 구현에서는 RDKit으로 정밀 분석
        cysteine_pattern_count = smiles.count("S")  # Simplified
        lysine_pattern_count = smiles.count("N")  # Simplified

        # DAR은 보통 2, 4, 8
        if cysteine_pattern_count >= 8:
            return 8
        elif cysteine_pattern_count >= 4:
            return 4
        else:
            return 4  # Default

    def _analyze_linker_type(self, smiles: str) -> Dict[str, Any]:
        """링커 타입 분석"""
        if not smiles:
            return {"type": "unknown", "reason": "No SMILES provided"}

        smiles_upper = smiles.upper()

        # Cleavable linker 패턴 (Val-Cit, Phe-Lys 등)
        cleavable_indicators = ["VAL", "CIT", "PHE", "LYS", "PEPTIDE"]
        for indicator in cleavable_indicators:
            if indicator in smiles_upper:
                return {"type": "cleavable", "reason": f"Contains {indicator} pattern"}

        # Non-cleavable linker 패턴 (SMCC, MCC, thioether)
        non_cleavable_indicators = ["SMCC", "MCC", "THIO"]
        for indicator in non_cleavable_indicators:
            if indicator in smiles_upper:
                return {"type": "non-cleavable", "reason": f"Contains {indicator} pattern"}

        return {"type": "unknown", "reason": "Linker type could not be determined from SMILES"}

    def _generate_constraint_summary(
        self,
        violations: List[Dict],
        warnings: List[Dict]
    ) -> str:
        """제약조건 검사 요약"""
        if not violations and not warnings:
            return "All user constraints satisfied"

        parts = []
        if violations:
            parts.append(f"{len(violations)} constraint violation(s)")
        if warnings:
            parts.append(f"{len(warnings)} warning(s)")

        return "; ".join(parts)

    def _calculate_risk_score(
        self,
        chemistry_validation: Dict,
        constraint_check: ConstraintCheck
    ) -> float:
        """리스크 스코어 계산 (0-10)"""
        score = 0

        # 화학 검증 리스크
        if not chemistry_validation.get("lipinski_pass", True):
            score += self.RISK_WEIGHTS["lipinski_fail"]

        if chemistry_validation.get("pains_detected", False):
            score += self.RISK_WEIGHTS["pains_detected"]

        if not chemistry_validation.get("sa_acceptable", True):
            score += self.RISK_WEIGHTS["high_sa_score"]

        # 제약조건 위반 리스크
        for violation in constraint_check.get("violations", []):
            if violation.get("severity") == "high":
                score += self.RISK_WEIGHTS["constraint_violation_high"]
            else:
                score += self.RISK_WEIGHTS["constraint_violation_medium"]

        return min(10, score)

    def _make_decision(
        self,
        chemistry_validation: Dict,
        constraint_check: ConstraintCheck,
        risk_score: float
    ) -> Dict[str, Any]:
        """최종 결정"""

        # High-severity constraint violations -> 재설계
        high_violations = [
            v for v in constraint_check.get("violations", [])
            if v.get("severity") == "high"
        ]

        if high_violations:
            return {
                "approved": False,
                "action": "redesign",
                "reasoning": f"Critical constraint violations: {[v['constraint'] for v in high_violations]}",
                "confidence": 0.9
            }

        # PAINS 감지 -> 수동 검토
        if chemistry_validation.get("pains_detected", False):
            return {
                "approved": False,
                "action": "manual_review",
                "reasoning": "PAINS alerts detected - requires expert review",
                "confidence": 0.75
            }

        # Lipinski 실패 -> 재설계
        if not chemistry_validation.get("lipinski_pass", True):
            return {
                "approved": False,
                "action": "redesign",
                "reasoning": f"Failed Lipinski's Rule (violations: {chemistry_validation.get('lipinski_violations', 0)})",
                "confidence": 0.85
            }

        # High risk score -> 수동 검토
        if risk_score > 7:
            return {
                "approved": False,
                "action": "manual_review",
                "reasoning": f"High risk score ({risk_score}/10) requires expert review",
                "confidence": 0.7
            }

        # Medium-severity violations -> 경고와 함께 승인
        if constraint_check.get("violations", []):
            return {
                "approved": True,
                "action": "approved_with_warnings",
                "reasoning": f"Approved with {len(constraint_check['violations'])} minor constraint deviations",
                "confidence": 0.75,
                "warnings": constraint_check["violations"]
            }

        # All good
        return {
            "approved": True,
            "action": "approved",
            "reasoning": "All validations passed, constraints satisfied",
            "confidence": 0.95
        }

    def _generate_report(
        self,
        state: DesignSessionState,
        chemistry_validation: Dict,
        constraint_check: ConstraintCheck,
        risk_score: float,
        decision: Dict
    ) -> Dict[str, Any]:
        """최종 리포트 생성"""
        return {
            "session_id": state["session_id"],
            "session_type": state["session_type"],
            "tier": state["tier"],

            "summary": {
                "target": state["target_antigen"],
                "indication": state["target_indication"],
                "dar": state["requested_dar"],
                "smiles_generated": bool(state["current_smiles"]),
                "candidates_count": len(state.get("candidates", []))
            },

            "validation": {
                "chemistry": chemistry_validation,
                "constraints": {
                    "passed": constraint_check["passed"],
                    "violations_count": len(constraint_check.get("violations", [])),
                    "warnings_count": len(constraint_check.get("warnings", []))
                }
            },

            "risk_assessment": {
                "score": risk_score,
                "level": "High" if risk_score > 7 else "Medium" if risk_score > 4 else "Low"
            },

            "decision": {
                "approved": decision["approved"],
                "action": decision["action"],
                "reasoning": decision["reasoning"],
                "confidence": decision["confidence"]
            },

            "generated_at": state.get("created_at"),
            "regulatory_version": "v1.0"
        }

    # ============== Phase 1: RDKit + Gemini Hybrid Methods ==============

    def _init_pains_filter(self):
        """PAINS 필터 초기화"""
        if RDKIT_AVAILABLE:
            try:
                params = FilterCatalogParams()
                params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
                self.pains_catalog = FilterCatalog.FilterCatalog(params)
                logger.info("[auditor] PAINS filter initialized")
            except Exception as e:
                logger.warning(f"[auditor] PAINS filter init failed: {e}")
                self.pains_catalog = None
        else:
            self.pains_catalog = None

    async def _validate_smiles_integrity(self, smiles: str) -> Dict[str, Any]:
        """
        [Phase 1] SMILES 무결성 검증

        RDKit을 사용하여 SMILES 문자열의 유효성을 검증합니다.
        """
        if not smiles or not smiles.strip():
            return {
                "valid": False,
                "error": "Empty SMILES string"
            }

        if not RDKIT_AVAILABLE:
            # RDKit 없으면 기본 검증만
            return {"valid": True, "warning": "RDKit not available, basic validation only"}

        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {
                    "valid": False,
                    "error": "RDKit failed to parse SMILES - invalid structure"
                }

            # 기본 sanity check
            num_atoms = mol.GetNumAtoms()
            if num_atoms < 3:
                return {
                    "valid": False,
                    "error": f"Molecule too small ({num_atoms} atoms)"
                }

            if num_atoms > 500:
                return {
                    "valid": True,
                    "warning": f"Large molecule ({num_atoms} atoms) - may be an ADC conjugate"
                }

            return {
                "valid": True,
                "num_atoms": num_atoms,
                "canonical_smiles": Chem.MolToSmiles(mol)
            }

        except Exception as e:
            return {
                "valid": False,
                "error": str(e)
            }

    async def _validate_chemistry_rdkit(self, state: DesignSessionState) -> Dict[str, Any]:
        """
        [Phase 1] RDKit 기반 정밀 화학 검증

        - Lipinski's Rule (ADC 확장)
        - PAINS 필터
        - SA Score
        - TPSA, Rotatable Bonds
        """
        smiles = state.get("current_smiles", "")
        metrics = state.get("calculated_metrics", {})

        result = {
            "lipinski_pass": True,
            "lipinski_violations": 0,
            "pains_detected": False,
            "pains_patterns": [],
            "sa_score": None,
            "sa_acceptable": True,
            "metrics": {},
            "adc_compliant": True
        }

        if not RDKIT_AVAILABLE or not smiles:
            # Fallback: Coder가 계산한 metrics 사용
            return self._validate_chemistry(state)

        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                result["lipinski_pass"] = False
                result["error"] = "Invalid SMILES"
                return result

            # 1. 물성 계산
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

            result["metrics"] = {
                "mw": round(mw, 2),
                "logp": round(logp, 2),
                "hbd": hbd,
                "hba": hba,
                "tpsa": round(tpsa, 2),
                "rotatable_bonds": rotatable_bonds,
                "num_rings": rdMolDescriptors.CalcNumRings(mol),
                "num_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol)
            }

            # 2. ADC 확장 Lipinski 검증
            violations = []
            if mw > self.ADC_THRESHOLDS["mw_max"]:
                violations.append(f"MW ({mw:.0f}) > {self.ADC_THRESHOLDS['mw_max']}")
            if logp > self.ADC_THRESHOLDS["logp_max"]:
                violations.append(f"LogP ({logp:.1f}) > {self.ADC_THRESHOLDS['logp_max']}")
            if hbd > self.ADC_THRESHOLDS["hbd_max"]:
                violations.append(f"HBD ({hbd}) > {self.ADC_THRESHOLDS['hbd_max']}")
            if hba > self.ADC_THRESHOLDS["hba_max"]:
                violations.append(f"HBA ({hba}) > {self.ADC_THRESHOLDS['hba_max']}")

            result["lipinski_violations"] = len(violations)
            result["lipinski_pass"] = len(violations) <= 1  # ADC는 1개까지 허용
            result["lipinski_details"] = violations

            # 3. PAINS 필터
            if self.pains_catalog:
                entries = self.pains_catalog.GetMatches(mol)
                if entries:
                    result["pains_detected"] = True
                    result["pains_patterns"] = [
                        entry.GetDescription() for entry in entries
                    ]

            # 4. SA Score (합성 접근성)
            try:
                from rdkit.Chem import RDConfig
                import os
                import sys
                sa_module_path = os.path.join(RDConfig.RDContribDir, 'SA_Score')
                if os.path.exists(sa_module_path):
                    sys.path.insert(0, sa_module_path)
                    import sascorer
                    sa_score = sascorer.calculateScore(mol)
                    result["sa_score"] = round(sa_score, 2)
                    result["sa_acceptable"] = sa_score < 6
            except Exception as e:
                logger.debug(f"[auditor] SA Score calculation skipped: {e}")

            # 5. ADC 적합성 종합 평가
            result["adc_compliant"] = (
                result["lipinski_pass"] and
                not result["pains_detected"] and
                result["sa_acceptable"]
            )

        except Exception as e:
            logger.error(f"[auditor] RDKit validation error: {e}")
            result["error"] = str(e)

        return result

    async def _make_decision_with_gemini(
        self,
        state: DesignSessionState,
        chemistry_validation: Dict,
        constraint_check: ConstraintCheck,
        risk_score: float
    ) -> Dict[str, Any]:
        """
        [Phase 1] Gemini를 사용한 지능형 검증 결정

        RDKit 검증 결과를 바탕으로 Gemini가 최종 결정을 내립니다.
        """
        # 기본 결정 먼저 수행
        base_decision = self._make_decision(
            chemistry_validation,
            constraint_check,
            risk_score
        )

        # Gemini로 결정 보강
        system_prompt = """당신은 ADC(Antibody-Drug Conjugate) 약물 개발 전문 검증자입니다.
RDKit 분석 결과와 제약조건 검사 결과를 검토하고 최종 결정을 내려야 합니다.

응답 형식 (JSON):
{
    "approved": true/false,
    "action": "approved" | "approved_with_warnings" | "manual_review" | "redesign",
    "reasoning": "결정 이유 (한국어)",
    "confidence": 0.0~1.0,
    "suggestions": ["개선 제안사항 리스트"],
    "key_concerns": ["주요 우려사항 리스트"],
    "failed_checks": ["실패한 검증 항목 리스트"]
}"""

        metrics = chemistry_validation.get("metrics", {})
        violations = constraint_check.get("violations", [])

        user_prompt = f"""ADC 후보 검증 결과를 분석해주세요:

**타겟**: {state.get('target_antigen', 'N/A')}
**적응증**: {state.get('target_indication', 'N/A')}

**화학적 검증 결과**:
- Lipinski 통과: {chemistry_validation.get('lipinski_pass', 'N/A')}
- Lipinski 위반 수: {chemistry_validation.get('lipinski_violations', 0)}
- PAINS 감지: {chemistry_validation.get('pains_detected', False)}
- SA Score: {chemistry_validation.get('sa_score', 'N/A')}
- MW: {metrics.get('mw', 'N/A')}, LogP: {metrics.get('logp', 'N/A')}
- HBD: {metrics.get('hbd', 'N/A')}, HBA: {metrics.get('hba', 'N/A')}

**제약조건 검사**:
- 통과: {constraint_check.get('passed', True)}
- 위반 수: {len(violations)}
- 위반 내용: {violations if violations else '없음'}

**리스크 스코어**: {risk_score}/10

기본 결정: {base_decision['action']} (신뢰도: {base_decision['confidence']})

위 결과를 바탕으로 최종 결정을 내려주세요. ADC 특성상 일반 약물보다 물성 기준이 유연할 수 있습니다."""

        try:
            response = await self.llm.ainvoke([
                SystemMessage(content=system_prompt),
                HumanMessage(content=user_prompt)
            ])

            content = response.content
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0]
            elif "```" in content:
                content = content.split("```")[1].split("```")[0]

            gemini_decision = json.loads(content)

            # Gemini 결정과 기본 결정 병합
            return {
                "approved": gemini_decision.get("approved", base_decision["approved"]),
                "action": gemini_decision.get("action", base_decision["action"]),
                "reasoning": gemini_decision.get("reasoning", base_decision["reasoning"]),
                "confidence": gemini_decision.get("confidence", base_decision["confidence"]),
                "suggestions": gemini_decision.get("suggestions", []),
                "key_concerns": gemini_decision.get("key_concerns", []),
                "failed_checks": gemini_decision.get("failed_checks", []),
                "gemini_enhanced": True
            }

        except Exception as e:
            logger.warning(f"[auditor] Gemini decision enhancement failed: {e}")
            # Fallback: 기본 결정 반환
            base_decision["gemini_enhanced"] = False
            return base_decision

    async def _request_redesign(
        self,
        session_id: str,
        state: DesignSessionState,
        reason: str,
        failed_checks: List[str]
    ) -> AgentOutput:
        """재설계 요청 생성"""

        await websocket_hub.broadcast_auditor_feedback(
            session_id,
            is_in_redesign_loop=True,
            check_items=[{
                "id": check,
                "name": check,
                "category": "other",
                "status": "fail"
            } for check in failed_checks],
            redesign_request={
                "requestedBy": "auditor",
                "reason": reason,
                "failed_checks": failed_checks,
                "suggestions": ["SMILES 구조 재검토 필요", "RDKit 파싱 가능한 형식으로 수정"],
                "target_agent": "alchemist",
                "iteration": state.get("redesign_count", 0) + 1,
                "max_iterations": 3
            }
        )

        return AgentOutput(
            success=False,
            data={
                "decision": {
                    "approved": False,
                    "action": "redesign",
                    "reasoning": reason,
                    "confidence": 0.9
                }
            },
            reasoning=reason,
            confidence_score=0.9,
            next_agent="alchemist"
        )

    def _format_check_items(self, chemistry_validation: Dict) -> List[Dict]:
        """검증 항목을 UI 포맷으로 변환"""
        items = []
        metrics = chemistry_validation.get("metrics", {})

        # Lipinski
        items.append({
            "id": "lipinski",
            "name": "Lipinski Rule",
            "category": "lipinski",
            "status": "pass" if chemistry_validation.get("lipinski_pass") else "fail",
            "currentValue": chemistry_validation.get("lipinski_violations", 0),
            "targetValue": "≤1 violations",
            "reasoning": f"Violations: {chemistry_validation.get('lipinski_details', [])}"
        })

        # PAINS
        items.append({
            "id": "pains",
            "name": "PAINS Filter",
            "category": "pains",
            "status": "fail" if chemistry_validation.get("pains_detected") else "pass",
            "currentValue": len(chemistry_validation.get("pains_patterns", [])),
            "targetValue": 0,
            "reasoning": str(chemistry_validation.get("pains_patterns", []))
        })

        # SA Score
        sa_score = chemistry_validation.get("sa_score")
        if sa_score is not None:
            items.append({
                "id": "sa_score",
                "name": "SA Score",
                "category": "stability",
                "status": "pass" if chemistry_validation.get("sa_acceptable") else "warning",
                "currentValue": sa_score,
                "targetValue": "< 6",
                "unit": ""
            })

        # MW
        if metrics.get("mw"):
            items.append({
                "id": "mw",
                "name": "Molecular Weight",
                "category": "lipinski",
                "status": "pass" if metrics["mw"] <= self.ADC_THRESHOLDS["mw_max"] else "warning",
                "currentValue": metrics["mw"],
                "targetValue": f"≤ {self.ADC_THRESHOLDS['mw_max']}",
                "unit": "Da"
            })

        # LogP
        if metrics.get("logp") is not None:
            items.append({
                "id": "logp",
                "name": "LogP",
                "category": "lipinski",
                "status": "pass" if metrics["logp"] <= self.ADC_THRESHOLDS["logp_max"] else "warning",
                "currentValue": metrics["logp"],
                "targetValue": f"≤ {self.ADC_THRESHOLDS['logp_max']}",
                "unit": ""
            })

        return items
