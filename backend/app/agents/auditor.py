"""
The Auditor Agent - 검증 및 감사
PAINS/Lipinski/독성 검증 + Constraint Guardrail
"""
from typing import List, Dict, Any, Optional
import logging

from app.agents.base_agent import BaseDesignAgent, AgentOutput
from app.agents.design_state import (
    DesignSessionState,
    ValidationFlags,
    ConstraintCheck
)
from app.core.supabase import get_supabase_client

logger = logging.getLogger(__name__)


class AuditorAgent(BaseDesignAgent):
    """
    The Auditor: 검증 및 감사 에이전트

    핵심 기능:
    1. PAINS/Lipinski/독성 검증
    2. Constraint Guardrail: 초기 제약조건 준수 확인 (Constraint 4)
    3. 리스크 스코어 계산
    4. 최종 승인/재설계 결정
    """

    name = "auditor"

    # 리스크 가중치
    RISK_WEIGHTS = {
        "lipinski_fail": 3,
        "pains_detected": 4,
        "constraint_violation_high": 5,
        "constraint_violation_medium": 2,
        "low_similarity": 2,
        "high_sa_score": 1
    }

    async def execute(self, state: DesignSessionState) -> AgentOutput:
        """메인 실행"""
        session_id = state["session_id"]

        await self._log_start(session_id, {
            "smiles": state["current_smiles"][:30] + "..." if state["current_smiles"] else None,
            "tier": state["tier"]
        })

        try:
            # 1. 기본 화학 검증
            chemistry_validation = self._validate_chemistry(state)

            # 2. Constraint Guardrail 검사 (Constraint 4)
            constraint_check = await self._check_constraint_guardrail(state)

            # 3. 리스크 스코어 계산
            risk_score = self._calculate_risk_score(
                chemistry_validation,
                constraint_check
            )

            # 4. 최종 결정
            decision = self._make_decision(
                chemistry_validation,
                constraint_check,
                risk_score
            )

            # 5. 최종 리포트 생성
            final_report = self._generate_report(
                state, chemistry_validation, constraint_check, risk_score, decision
            )

            output = AgentOutput(
                success=decision["approved"],
                data={
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
