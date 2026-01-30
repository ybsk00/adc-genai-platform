"""
NVIDIA NIM Service with Adaptive API Interface
Phase 1: Gemini 3.0 Pro Fallback
Phase 3: NVIDIA BioNeMo NIM Production

어댑터 패턴을 적용하여 USE_NIM_API 플래그에 따라 자동 스위칭
"""
import json
import platform
import hashlib
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List, Callable
from datetime import datetime
import httpx

from app.core.config import settings
from app.core.gemini import get_gemini_model


class StructurePredictorInterface(ABC):
    """구조 예측 인터페이스 (어댑터 패턴)"""

    @abstractmethod
    async def predict_structure(
        self,
        sequence: str,
        ligand_smiles: Optional[str] = None,
        progress_callback: Optional[Callable] = None
    ) -> Dict[str, Any]:
        """항체/단백질 구조 예측"""
        pass

    @abstractmethod
    async def run_msa(self, sequence: str) -> Dict[str, Any]:
        """Multiple Sequence Alignment 생성"""
        pass

    @abstractmethod
    async def generate_molecule(
        self,
        scaffold: str,
        constraints: Dict[str, Any]
    ) -> Dict[str, Any]:
        """분자 생성 (MolMIM)"""
        pass

    @abstractmethod
    async def predict_toxicity(
        self,
        smiles: str,
        target: str
    ) -> Dict[str, Any]:
        """독성 예측"""
        pass


class GeminiFallbackPredictor(StructurePredictorInterface):
    """
    Phase 1: Gemini 3.0 Pro 기반 Fallback 예측기

    NVIDIA NIM 연동 전, 모든 에이전트가 외부 API 장애 없이
    로직을 테스트할 수 있도록 Gemini의 SOTA 추론 능력 활용
    """

    def __init__(self):
        self.model = get_gemini_model(
            temperature=0.1,
            model_name=settings.GEMINI_PRO_MODEL_ID
        )
        self.mode = "fallback"
        self.engine_name = "Gemini 3.0 Pro Preview"
        self.engine_icon = "gemini-3-pro-preview"

    async def predict_structure(
        self,
        sequence: str,
        ligand_smiles: Optional[str] = None,
        progress_callback: Optional[Callable] = None
    ) -> Dict[str, Any]:
        """
        [Fallback] Gemini 3 Pro Preview 활용
        구조 예측 대신 'Vibe-based' 논리적 구조 분석 및 Placeholder 반환
        """

        if progress_callback:
            await progress_callback(10, "Gemini 추론 엔진 초기화 중...")

        prompt = f"""You are an expert structural biologist and antibody engineer.
Analyze the following antibody/protein sequence and provide a logical structural prediction
based on known antibody architectures and structural biology principles.

Sequence ({len(sequence)} residues):
{sequence[:500]}{'...' if len(sequence) > 500 else ''}

{f'Ligand SMILES: {ligand_smiles}' if ligand_smiles else ''}

Provide a comprehensive structural analysis in JSON format:
{{
    "domain_structure": {{
        "domains": ["VH", "VL", "CH1", "CH2", "CH3", "CL"],
        "boundaries": {{"VH": [1, 120], "VL": [121, 230], ...}},
        "confidence": 0.85
    }},
    "cdr_regions": {{
        "CDR-H1": {{"start": 26, "end": 35, "sequence": "..."}},
        "CDR-H2": {{"start": 50, "end": 65, "sequence": "..."}},
        "CDR-H3": {{"start": 93, "end": 102, "sequence": "..."}},
        "CDR-L1": {{"start": 24, "end": 34, "sequence": "..."}},
        "CDR-L2": {{"start": 50, "end": 56, "sequence": "..."}},
        "CDR-L3": {{"start": 89, "end": 97, "sequence": "..."}}
    }},
    "conjugation_sites": {{
        "cysteine_positions": [239, 242, 280, 318],
        "lysine_positions": [246, 288, 290, 317, 320],
        "recommended_sites": ["Cys-239", "Cys-242"],
        "dar_potential": 4
    }},
    "stability_prediction": {{
        "overall_score": 85,
        "aggregation_risk": "Low",
        "thermal_stability": "High",
        "ph_sensitivity": "Moderate"
    }},
    "structural_features": {{
        "glycosylation_sites": ["N297"],
        "disulfide_bonds": ["C22-C96", "C148-C204", "C226-C229"],
        "unusual_features": []
    }},
    "rationale": "Detailed explanation of the structural analysis..."
}}

Return ONLY valid JSON, no additional text."""

        if progress_callback:
            await progress_callback(30, "구조 분석 중...")

        try:
            response = await self.model.ainvoke(prompt)
            result_text = response.content

            # JSON 파싱 시도
            try:
                # JSON 블록 추출
                if "```json" in result_text:
                    result_text = result_text.split("```json")[1].split("```")[0]
                elif "```" in result_text:
                    result_text = result_text.split("```")[1].split("```")[0]

                analysis = json.loads(result_text.strip())
            except json.JSONDecodeError:
                analysis = {"raw_response": result_text, "parse_error": True}

            if progress_callback:
                await progress_callback(90, "결과 정리 중...")

            return {
                "source": "gemini-3-pro-preview",
                "mode": "fallback",
                "prediction_type": "logical_inference",
                "sequence_length": len(sequence),
                "analysis": analysis,
                "domain_structure": analysis.get("domain_structure", {}),
                "cdr_regions": analysis.get("cdr_regions", {}),
                "conjugation_sites": analysis.get("conjugation_sites", {}).get("recommended_sites", []),
                "stability_score": analysis.get("stability_prediction", {}).get("overall_score", 75),
                "aggregation_risk": analysis.get("stability_prediction", {}).get("aggregation_risk", "Medium"),
                "pdb_data": None,  # Placeholder - Phase 3에서 실제 PDB
                "plddt_mean": None,  # Placeholder
                "confidence": 0.75,
                "engine_icon": self.engine_icon,
                "_warning": "Preview Mode: 실제 3D 구조는 NVIDIA NIM 연동 후 제공됩니다.",
                "_preview_mode": True
            }

        except Exception as e:
            return {
                "source": "gemini-3-pro-preview",
                "mode": "fallback",
                "error": str(e),
                "confidence": 0,
                "engine_icon": self.engine_icon
            }

    async def run_msa(self, sequence: str) -> Dict[str, Any]:
        """
        [Fallback] MSA 생성 대신 논리적 시퀀스 분석
        """

        prompt = f"""Analyze this antibody/protein sequence for evolutionary and structural context.

Sequence ({len(sequence)} residues):
{sequence[:300]}{'...' if len(sequence) > 300 else ''}

Provide analysis in JSON format:
{{
    "sequence_type": "IgG1/IgG2/IgG4/other",
    "species_origin": "human/humanized/chimeric/murine",
    "framework_regions": {{
        "FR1": {{"conservation": "high/medium/low"}},
        "FR2": {{"conservation": "high/medium/low"}},
        "FR3": {{"conservation": "high/medium/low"}},
        "FR4": {{"conservation": "high/medium/low"}}
    }},
    "variable_regions": {{
        "vh_family": "VH1/VH3/VH4/other",
        "vl_family": "Vkappa1/Vlambda1/other"
    }},
    "evolutionary_context": {{
        "closest_germline": "IGHV3-23",
        "somatic_mutations": 15,
        "affinity_maturation_level": "high/medium/low"
    }},
    "therapeutic_potential": {{
        "immunogenicity_risk": "low/medium/high",
        "developability_score": 85
    }},
    "rationale": "Explanation of the analysis..."
}}

Return ONLY valid JSON."""

        try:
            response = await self.model.ainvoke(prompt)
            result_text = response.content

            try:
                if "```json" in result_text:
                    result_text = result_text.split("```json")[1].split("```")[0]
                analysis = json.loads(result_text.strip())
            except json.JSONDecodeError:
                analysis = {"raw_response": result_text}

            return {
                "source": "gemini-3-pro-preview",
                "mode": "fallback",
                "num_sequences": 0,  # Placeholder - 실제 MSA 미생성
                "analysis": analysis,
                "sequence_type": analysis.get("sequence_type", "Unknown"),
                "framework_conservation": analysis.get("framework_regions", {}),
                "evolutionary_context": analysis.get("evolutionary_context", {}),
                "engine_icon": self.engine_icon,
                "_warning": "Preview Mode: 실제 MSA는 NVIDIA NIM 연동 후 제공됩니다.",
                "_preview_mode": True
            }

        except Exception as e:
            return {
                "source": "gemini-3-pro-preview",
                "mode": "fallback",
                "error": str(e),
                "engine_icon": self.engine_icon
            }

    async def generate_molecule(
        self,
        scaffold: str,
        constraints: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        [Fallback] 분자 생성 대신 Gemini 기반 논리적 설계
        """

        prompt = f"""You are an expert medicinal chemist specializing in ADC (Antibody-Drug Conjugate) design.
Based on the following scaffold and constraints, suggest optimized molecular structures.

Scaffold SMILES: {scaffold}

Constraints:
- Target: {constraints.get('target', 'Unknown')}
- DAR (Drug-to-Antibody Ratio): {constraints.get('dar', 4)}
- Linker Type: {constraints.get('linker_type', 'cleavable')}
- Design Goal: {constraints.get('design_goal', 'Optimize efficacy and safety')}

Generate 3-5 candidate molecules in JSON format:
{{
    "candidates": [
        {{
            "smiles": "valid SMILES string",
            "name": "Candidate-1",
            "modifications": ["description of modifications from scaffold"],
            "predicted_properties": {{
                "mw": 450,
                "logp": 2.5,
                "hbd": 2,
                "hba": 5,
                "psa": 80
            }},
            "rationale": "Why this modification improves the design",
            "score": 0.85
        }}
    ],
    "design_strategy": "Overall design rationale",
    "warnings": ["Any concerns about the designs"]
}}

Return ONLY valid JSON with chemically valid SMILES."""

        try:
            response = await self.model.ainvoke(prompt)
            result_text = response.content

            try:
                if "```json" in result_text:
                    result_text = result_text.split("```json")[1].split("```")[0]
                result = json.loads(result_text.strip())
            except json.JSONDecodeError:
                result = {"candidates": [], "raw_response": result_text}

            candidates = result.get("candidates", [])

            return {
                "source": "gemini-3-pro-preview",
                "mode": "fallback",
                "candidates": candidates,
                "design_strategy": result.get("design_strategy", ""),
                "warnings": result.get("warnings", []),
                "num_candidates": len(candidates),
                "engine_icon": self.engine_icon,
                "_warning": "Preview Mode: 실제 MolMIM 생성은 NVIDIA NIM 연동 후 제공됩니다.",
                "_preview_mode": True
            }

        except Exception as e:
            return {
                "source": "gemini-3-pro-preview",
                "mode": "fallback",
                "error": str(e),
                "candidates": [],
                "engine_icon": self.engine_icon
            }

    async def predict_toxicity(
        self,
        smiles: str,
        target: str
    ) -> Dict[str, Any]:
        """
        [Fallback] 독성 예측
        """

        prompt = f"""You are a toxicologist specializing in drug safety assessment.
Analyze the following molecule for potential toxicity concerns in the context of ADC therapy.

SMILES: {smiles}
Target: {target}

Provide toxicity analysis in JSON format:
{{
    "overall_risk_score": 0.35,
    "risk_level": "Low/Medium/High",
    "toxicity_endpoints": {{
        "hepatotoxicity": {{"risk": "low", "confidence": 0.8}},
        "cardiotoxicity": {{"risk": "low", "confidence": 0.75}},
        "nephrotoxicity": {{"risk": "medium", "confidence": 0.7}},
        "neurotoxicity": {{"risk": "low", "confidence": 0.85}},
        "hematotoxicity": {{"risk": "medium", "confidence": 0.8}},
        "immunogenicity": {{"risk": "low", "confidence": 0.7}}
    }},
    "structural_alerts": [
        {{"alert": "description", "severity": "low/medium/high"}}
    ],
    "off_target_concerns": [
        {{"target": "protein name", "binding_probability": 0.2}}
    ],
    "recommendations": ["safety recommendations"],
    "rationale": "Explanation of the toxicity assessment"
}}

Return ONLY valid JSON."""

        try:
            response = await self.model.ainvoke(prompt)
            result_text = response.content

            try:
                if "```json" in result_text:
                    result_text = result_text.split("```json")[1].split("```")[0]
                analysis = json.loads(result_text.strip())
            except json.JSONDecodeError:
                analysis = {"raw_response": result_text}

            return {
                "source": "gemini-3-pro-preview",
                "mode": "fallback",
                "smiles": smiles,
                "target": target,
                "overall_risk_score": analysis.get("overall_risk_score", 0.5),
                "risk_level": analysis.get("risk_level", "Medium"),
                "toxicity_endpoints": analysis.get("toxicity_endpoints", {}),
                "structural_alerts": analysis.get("structural_alerts", []),
                "off_target_concerns": analysis.get("off_target_concerns", []),
                "recommendations": analysis.get("recommendations", []),
                "engine_icon": self.engine_icon,
                "_warning": "Preview Mode: 실제 구조 기반 독성 분석은 NVIDIA NIM 연동 후 제공됩니다.",
                "_preview_mode": True
            }

        except Exception as e:
            return {
                "source": "gemini-3-pro-preview",
                "mode": "fallback",
                "error": str(e),
                "engine_icon": self.engine_icon
            }


class NVIDIANIMPredictor(StructurePredictorInterface):
    """
    Phase 3: 실제 NVIDIA BioNeMo NIM API
    AlphaFold 3, MolMIM, ESMFold 등 연동
    """

    def __init__(self, api_key: str):
        self.api_key = api_key
        self.base_url = settings.NVIDIA_NIM_BASE_URL
        self.mode = "production"
        self.engine_name = "NVIDIA BioNeMo NIM"
        self.engine_icon = "nvidia-bionemo"
        self.client = httpx.AsyncClient(timeout=300)  # 5분 타임아웃

    async def _call_api(
        self,
        endpoint: str,
        payload: Dict[str, Any]
    ) -> Dict[str, Any]:
        """NIM API 호출"""

        response = await self.client.post(
            f"{self.base_url}/{endpoint}",
            headers={
                "Authorization": f"Bearer {self.api_key}",
                "Content-Type": "application/json"
            },
            json=payload
        )

        response.raise_for_status()
        return response.json()

    async def predict_structure(
        self,
        sequence: str,
        ligand_smiles: Optional[str] = None,
        progress_callback: Optional[Callable] = None
    ) -> Dict[str, Any]:
        """실제 AlphaFold 3 구조 예측"""

        if progress_callback:
            await progress_callback(10, "NVIDIA NIM 서버 연결 중...")

        try:
            payload = {
                "sequences": [sequence],
                "ligands": [{"smiles": ligand_smiles}] if ligand_smiles else []
            }

            if progress_callback:
                await progress_callback(30, "AlphaFold 3 구조 예측 중...")

            result = await self._call_api("alphafold3/predict", payload)

            if progress_callback:
                await progress_callback(90, "결과 처리 중...")

            return {
                "source": "nvidia-bionemo-nim",
                "mode": "production",
                "prediction_type": "alphafold3",
                "sequence_length": len(sequence),
                "pdb_data": result.get("pdb_string"),
                "plddt_mean": result.get("plddt_mean"),
                "confidence": result.get("confidence", 0.9),
                "conjugation_sites": result.get("conjugation_sites", []),
                "stability_score": result.get("stability_score"),
                "engine_icon": self.engine_icon,
                "_preview_mode": False
            }

        except Exception as e:
            return {
                "source": "nvidia-bionemo-nim",
                "mode": "production",
                "error": str(e),
                "engine_icon": self.engine_icon
            }

    async def run_msa(self, sequence: str) -> Dict[str, Any]:
        """실제 고속 MSA 생성 (40분 → 1.5분)"""

        try:
            result = await self._call_api("esmfold/msa", {
                "sequence": sequence
            })

            return {
                "source": "nvidia-bionemo-nim",
                "mode": "production",
                "num_sequences": result.get("num_sequences", 0),
                "msa_data": result.get("msa"),
                "processing_time_seconds": result.get("processing_time"),
                "engine_icon": self.engine_icon,
                "_preview_mode": False
            }

        except Exception as e:
            return {
                "source": "nvidia-bionemo-nim",
                "mode": "production",
                "error": str(e),
                "engine_icon": self.engine_icon
            }

    async def generate_molecule(
        self,
        scaffold: str,
        constraints: Dict[str, Any]
    ) -> Dict[str, Any]:
        """실제 MolMIM 분자 생성"""

        try:
            result = await self._call_api("molmim/generate", {
                "scaffold": scaffold,
                "constraints": constraints,
                "num_candidates": constraints.get("num_candidates", 5)
            })

            return {
                "source": "nvidia-bionemo-nim",
                "mode": "production",
                "candidates": result.get("candidates", []),
                "generation_time_seconds": result.get("processing_time"),
                "engine_icon": self.engine_icon,
                "_preview_mode": False
            }

        except Exception as e:
            return {
                "source": "nvidia-bionemo-nim",
                "mode": "production",
                "error": str(e),
                "candidates": [],
                "engine_icon": self.engine_icon
            }

    async def predict_toxicity(
        self,
        smiles: str,
        target: str
    ) -> Dict[str, Any]:
        """구조 기반 독성 예측"""

        try:
            result = await self._call_api("toxicity/predict", {
                "smiles": smiles,
                "target": target
            })

            return {
                "source": "nvidia-bionemo-nim",
                "mode": "production",
                "overall_risk_score": result.get("risk_score"),
                "risk_level": result.get("risk_level"),
                "toxicity_endpoints": result.get("endpoints", {}),
                "engine_icon": self.engine_icon,
                "_preview_mode": False
            }

        except Exception as e:
            return {
                "source": "nvidia-bionemo-nim",
                "mode": "production",
                "error": str(e),
                "engine_icon": self.engine_icon
            }


class NVIDIANIMService:
    """
    Adaptive API Interface

    환경변수 USE_NIM_API에 따라 Gemini Fallback 또는 실제 NIM API 사용
    Phase 1: Gemini 3.0 Pro Fallback (디버깅 모드)
    Phase 3: NVIDIA BioNeMo NIM (프로덕션 모드)
    """

    def __init__(self):
        # 환경변수에 따른 스위칭 로직
        self.use_nim_api = settings.USE_NIM_API

        if self.use_nim_api and settings.NVIDIA_NIM_API_KEY:
            self.predictor = NVIDIANIMPredictor(settings.NVIDIA_NIM_API_KEY)
            self.engine_mode = "nvidia-nim"
            self.engine_icon = "nvidia-bionemo"
            self.engine_name = "NVIDIA BioNeMo NIM"
            self.status_message = "NVIDIA BioNeMo NIM 연동 중"
        else:
            self.predictor = GeminiFallbackPredictor()
            self.engine_mode = "gemini-fallback"
            self.engine_icon = "gemini-3-pro-preview"
            self.engine_name = "Gemini 3.0 Pro Preview"
            self.status_message = "SOTA 추론 엔진 가동 중"

    async def predict_structure(
        self,
        sequence: str,
        ligand_smiles: Optional[str] = None,
        progress_callback: Optional[Callable] = None
    ) -> Dict[str, Any]:
        """구조 예측 (어댑터 자동 스위칭)"""

        result = await self.predictor.predict_structure(
            sequence=sequence,
            ligand_smiles=ligand_smiles,
            progress_callback=progress_callback
        )
        result["engine_mode"] = self.engine_mode
        result["engine_name"] = self.engine_name
        return result

    async def run_msa(self, sequence: str) -> Dict[str, Any]:
        """MSA 생성 (어댑터 자동 스위칭)"""

        result = await self.predictor.run_msa(sequence)
        result["engine_mode"] = self.engine_mode
        result["engine_name"] = self.engine_name
        return result

    async def generate_molecule(
        self,
        scaffold: str,
        constraints: Dict[str, Any]
    ) -> Dict[str, Any]:
        """분자 생성 (어댑터 자동 스위칭)"""

        result = await self.predictor.generate_molecule(scaffold, constraints)
        result["engine_mode"] = self.engine_mode
        result["engine_name"] = self.engine_name
        return result

    async def predict_toxicity(
        self,
        smiles: str,
        target: str
    ) -> Dict[str, Any]:
        """독성 예측 (어댑터 자동 스위칭)"""

        result = await self.predictor.predict_toxicity(smiles, target)
        result["engine_mode"] = self.engine_mode
        result["engine_name"] = self.engine_name
        return result

    def get_engine_status(self) -> Dict[str, Any]:
        """
        현재 엔진 상태 반환 (UI 표시용)

        Returns:
            엔진 상태 정보 (mode, icon, status_message 등)
        """
        return {
            "mode": self.engine_mode,
            "icon": self.engine_icon,
            "name": self.engine_name,
            "use_nim_api": self.use_nim_api,
            "status_message": self.status_message,
            "is_preview_mode": not self.use_nim_api,
            "capabilities": {
                "structure_prediction": True,
                "msa_generation": True,
                "molecule_generation": True,
                "toxicity_prediction": True
            }
        }

    def get_node_id(self) -> Optional[str]:
        """NVIDIA NIM 노드 ID 반환 (Digital Lineage용)"""
        if self.use_nim_api:
            return "nim-production-node"  # 실제로는 API에서 반환
        return None


# Singleton instance
nvidia_nim_service = NVIDIANIMService()
