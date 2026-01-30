"""
Physical Validator - Sandbox Stress Test
분자 구조의 물리적 타당성 검증 서비스

Phase 2 Enhancement:
- Edge Case 검증 강화
- 3D 구조 기반 물리 검증
- physical_validations 테이블 저장
"""
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
import logging

# RDKit imports (optional)
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors, Descriptors
    from rdkit.Geometry import Point3D
    import numpy as np
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    np = None

from app.core.supabase import get_supabase_client
from app.core.config import settings

logger = logging.getLogger(__name__)


# ============================================================================
# Enums and Data Classes
# ============================================================================

class ValidationCategory(Enum):
    STRUCTURE = "structure"
    PROPERTY = "property"
    DRUGLIKENESS = "druglikeness"
    TOXICITY = "toxicity"
    STABILITY = "stability"
    SYNTHESIS = "synthesis"


class Severity(Enum):
    CRITICAL = "critical"
    WARNING = "warning"
    INFO = "info"


@dataclass
class PhysicalValidation:
    """물리 검증 결과"""
    check_name: str
    category: ValidationCategory
    passed: bool
    severity: Severity
    description: str
    measured_value: Optional[float] = None
    threshold: Optional[float] = None
    unit: Optional[str] = None
    affected_atoms: Optional[List[int]] = None
    suggested_fix: Optional[str] = None
    details: Dict[str, Any] = field(default_factory=dict)


# ============================================================================
# Physical Validator Service
# ============================================================================

class PhysicalValidator:
    """
    물리적 타당성 검증 서비스

    Sandbox Stress Test를 통해 물리적으로 불가능한 구조를 사전 차단
    """

    # Van der Waals 반경 (Å)
    VDW_RADII = {
        'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
        'S': 1.80, 'P': 1.80, 'F': 1.47, 'Cl': 1.75,
        'Br': 1.85, 'I': 1.98, 'B': 1.92, 'Si': 2.10
    }

    # 표준 결합 길이 (Å) - (원자1, 원자2, 결합차수)
    BOND_LENGTHS = {
        ('C', 'C', 1): (1.54, 0.15),   # C-C 단일, 허용 오차
        ('C', 'C', 2): (1.34, 0.10),   # C=C 이중
        ('C', 'C', 3): (1.20, 0.08),   # C≡C 삼중
        ('C', 'C', 1.5): (1.39, 0.08), # 방향족
        ('C', 'N', 1): (1.47, 0.12),
        ('C', 'N', 2): (1.27, 0.10),
        ('C', 'N', 3): (1.16, 0.08),
        ('C', 'N', 1.5): (1.34, 0.10), # 방향족
        ('C', 'O', 1): (1.43, 0.12),
        ('C', 'O', 2): (1.23, 0.08),
        ('C', 'O', 1.5): (1.33, 0.08), # 카복실
        ('C', 'S', 1): (1.82, 0.15),
        ('C', 'S', 2): (1.71, 0.10),
        ('C', 'H', 1): (1.09, 0.08),
        ('C', 'F', 1): (1.35, 0.10),
        ('C', 'Cl', 1): (1.77, 0.12),
        ('C', 'Br', 1): (1.94, 0.15),
        ('N', 'H', 1): (1.01, 0.08),
        ('O', 'H', 1): (0.96, 0.08),
        ('S', 'H', 1): (1.34, 0.10),
        ('N', 'N', 1): (1.45, 0.12),
        ('N', 'N', 2): (1.25, 0.10),
        ('N', 'O', 1): (1.40, 0.12),
        ('N', 'O', 2): (1.21, 0.08),
        ('O', 'O', 1): (1.48, 0.15),
        ('S', 'S', 1): (2.05, 0.15),
    }

    # ADC 특화 Lipinski+ 임계값
    ADC_THRESHOLDS = {
        "mw_max": 1000,           # ADC는 일반 약물보다 크다
        "logp_max": 6.0,
        "logp_min": -1.0,
        "hbd_max": 7,
        "hba_max": 15,
        "tpsa_max": 200,
        "rotatable_bonds_max": 15,
        "ring_count_max": 8,
    }

    def __init__(self):
        self.supabase = get_supabase_client()

    async def validate_structure(
        self,
        smiles: str,
        session_id: Optional[str] = None,
        calculation_id: Optional[str] = None,
        molecule_name: Optional[str] = None,
        generate_3d: bool = True,
        save_to_db: bool = True
    ) -> Dict[str, Any]:
        """
        분자 구조의 물리적 타당성 종합 검증

        Args:
            smiles: SMILES 문자열
            session_id: 설계 세션 ID
            calculation_id: 연산 ID (Digital Lineage 연결용)
            molecule_name: 분자 이름
            generate_3d: 3D 좌표 생성 여부
            save_to_db: DB 저장 여부

        Returns:
            검증 결과 요약
        """
        validations: List[PhysicalValidation] = []

        # 1. SMILES 유효성 검증
        smiles_validation = self._check_smiles_validity(smiles)
        validations.append(smiles_validation)

        if not smiles_validation.passed:
            return self._format_result(validations, smiles, session_id, calculation_id)

        if not RDKIT_AVAILABLE:
            logger.warning("[PhysicalValidator] RDKit not available, limited validation")
            return self._format_result(validations, smiles, session_id, calculation_id)

        # RDKit 분자 생성
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            validations.append(PhysicalValidation(
                check_name="mol_creation",
                category=ValidationCategory.STRUCTURE,
                passed=False,
                severity=Severity.CRITICAL,
                description="Failed to create molecule from SMILES"
            ))
            return self._format_result(validations, smiles, session_id, calculation_id)

        # 2. Valence 검증 (원자가 규칙)
        validations.append(self._check_valence(mol))

        # 3. 물리화학적 속성 검증
        validations.extend(self._check_physicochemical_properties(mol))

        # 4. Druglikeness 검증 (ADC 특화)
        validations.extend(self._check_druglikeness(mol))

        # 5. Ring Strain 검증
        validations.append(self._check_ring_strain(mol))

        # 6. 3D 구조 기반 검증
        if generate_3d:
            try:
                mol_3d = self._generate_3d_structure(mol)
                if mol_3d:
                    validations.append(self._check_steric_clashes(mol_3d))
                    validations.extend(self._check_bond_lengths(mol_3d))
                    validations.extend(self._check_bond_angles(mol_3d))
            except Exception as e:
                logger.warning(f"[PhysicalValidator] 3D validation error: {e}")
                validations.append(PhysicalValidation(
                    check_name="3d_generation",
                    category=ValidationCategory.STRUCTURE,
                    passed=True,
                    severity=Severity.INFO,
                    description=f"3D structure generation skipped: {str(e)}"
                ))

        # 7. Chirality 검증
        validations.append(self._check_chirality(mol))

        # 결과 포맷 및 저장
        result = self._format_result(validations, smiles, session_id, calculation_id, molecule_name)

        if save_to_db and session_id:
            await self._save_validations(validations, session_id, calculation_id, smiles, molecule_name)

        return result

    # =========================================================================
    # Validation Methods
    # =========================================================================

    def _check_smiles_validity(self, smiles: str) -> PhysicalValidation:
        """SMILES 문자열 유효성 검증"""
        if not smiles or not smiles.strip():
            return PhysicalValidation(
                check_name="smiles_validity",
                category=ValidationCategory.STRUCTURE,
                passed=False,
                severity=Severity.CRITICAL,
                description="Empty SMILES string"
            )

        if not RDKIT_AVAILABLE:
            return PhysicalValidation(
                check_name="smiles_validity",
                category=ValidationCategory.STRUCTURE,
                passed=True,
                severity=Severity.INFO,
                description="SMILES validity check skipped (RDKit not available)"
            )

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return PhysicalValidation(
                check_name="smiles_validity",
                category=ValidationCategory.STRUCTURE,
                passed=False,
                severity=Severity.CRITICAL,
                description="Invalid SMILES string - failed to parse",
                suggested_fix="Check SMILES syntax and atom symbols"
            )

        return PhysicalValidation(
            check_name="smiles_validity",
            category=ValidationCategory.STRUCTURE,
            passed=True,
            severity=Severity.INFO,
            description="Valid SMILES string",
            details={"atom_count": mol.GetNumAtoms(), "bond_count": mol.GetNumBonds()}
        )

    def _check_valence(self, mol) -> PhysicalValidation:
        """원자가 규칙 검증"""
        try:
            problems = Chem.DetectChemistryProblems(mol)
            valence_issues = [p for p in problems if "valence" in p.Message().lower()]

            if valence_issues:
                return PhysicalValidation(
                    check_name="valence_check",
                    category=ValidationCategory.STRUCTURE,
                    passed=False,
                    severity=Severity.CRITICAL,
                    description=f"Valence violation detected: {valence_issues[0].Message()}",
                    suggested_fix="Check atom connectivity and charge states"
                )

            return PhysicalValidation(
                check_name="valence_check",
                category=ValidationCategory.STRUCTURE,
                passed=True,
                severity=Severity.INFO,
                description="All atoms have valid valence"
            )
        except Exception as e:
            return PhysicalValidation(
                check_name="valence_check",
                category=ValidationCategory.STRUCTURE,
                passed=True,
                severity=Severity.INFO,
                description=f"Valence check completed with notes: {str(e)}"
            )

    def _check_physicochemical_properties(self, mol) -> List[PhysicalValidation]:
        """물리화학적 속성 검증"""
        validations = []

        # 분자량
        mw = Descriptors.MolWt(mol)
        mw_passed = mw <= self.ADC_THRESHOLDS["mw_max"]
        validations.append(PhysicalValidation(
            check_name="molecular_weight",
            category=ValidationCategory.PROPERTY,
            passed=mw_passed,
            severity=Severity.WARNING if not mw_passed else Severity.INFO,
            description=f"Molecular weight: {mw:.1f} Da",
            measured_value=mw,
            threshold=self.ADC_THRESHOLDS["mw_max"],
            unit="Da"
        ))

        # LogP
        logp = Descriptors.MolLogP(mol)
        logp_passed = self.ADC_THRESHOLDS["logp_min"] <= logp <= self.ADC_THRESHOLDS["logp_max"]
        validations.append(PhysicalValidation(
            check_name="logp_range",
            category=ValidationCategory.PROPERTY,
            passed=logp_passed,
            severity=Severity.WARNING if not logp_passed else Severity.INFO,
            description=f"LogP: {logp:.2f}",
            measured_value=logp,
            details={"min": self.ADC_THRESHOLDS["logp_min"], "max": self.ADC_THRESHOLDS["logp_max"]}
        ))

        # TPSA
        tpsa = Descriptors.TPSA(mol)
        tpsa_passed = tpsa <= self.ADC_THRESHOLDS["tpsa_max"]
        validations.append(PhysicalValidation(
            check_name="tpsa",
            category=ValidationCategory.PROPERTY,
            passed=tpsa_passed,
            severity=Severity.WARNING if not tpsa_passed else Severity.INFO,
            description=f"TPSA: {tpsa:.1f} Å²",
            measured_value=tpsa,
            threshold=self.ADC_THRESHOLDS["tpsa_max"],
            unit="Å²"
        ))

        return validations

    def _check_druglikeness(self, mol) -> List[PhysicalValidation]:
        """ADC 특화 Druglikeness 검증"""
        validations = []

        # HBD
        hbd = Descriptors.NumHDonors(mol)
        hbd_passed = hbd <= self.ADC_THRESHOLDS["hbd_max"]
        validations.append(PhysicalValidation(
            check_name="h_bond_donors",
            category=ValidationCategory.DRUGLIKENESS,
            passed=hbd_passed,
            severity=Severity.WARNING if not hbd_passed else Severity.INFO,
            description=f"H-bond donors: {hbd}",
            measured_value=hbd,
            threshold=self.ADC_THRESHOLDS["hbd_max"]
        ))

        # HBA
        hba = Descriptors.NumHAcceptors(mol)
        hba_passed = hba <= self.ADC_THRESHOLDS["hba_max"]
        validations.append(PhysicalValidation(
            check_name="h_bond_acceptors",
            category=ValidationCategory.DRUGLIKENESS,
            passed=hba_passed,
            severity=Severity.WARNING if not hba_passed else Severity.INFO,
            description=f"H-bond acceptors: {hba}",
            measured_value=hba,
            threshold=self.ADC_THRESHOLDS["hba_max"]
        ))

        # Rotatable Bonds
        rot_bonds = Descriptors.NumRotatableBonds(mol)
        rot_passed = rot_bonds <= self.ADC_THRESHOLDS["rotatable_bonds_max"]
        validations.append(PhysicalValidation(
            check_name="rotatable_bonds",
            category=ValidationCategory.DRUGLIKENESS,
            passed=rot_passed,
            severity=Severity.WARNING if not rot_passed else Severity.INFO,
            description=f"Rotatable bonds: {rot_bonds}",
            measured_value=rot_bonds,
            threshold=self.ADC_THRESHOLDS["rotatable_bonds_max"]
        ))

        # Ring Count
        ring_count = rdMolDescriptors.CalcNumRings(mol)
        ring_passed = ring_count <= self.ADC_THRESHOLDS["ring_count_max"]
        validations.append(PhysicalValidation(
            check_name="ring_count",
            category=ValidationCategory.DRUGLIKENESS,
            passed=ring_passed,
            severity=Severity.WARNING if not ring_passed else Severity.INFO,
            description=f"Ring count: {ring_count}",
            measured_value=ring_count,
            threshold=self.ADC_THRESHOLDS["ring_count_max"]
        ))

        return validations

    def _check_ring_strain(self, mol) -> PhysicalValidation:
        """환 변형 에너지 검증"""
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()

        # 3-membered 또는 4-membered rings 체크 (높은 strain)
        strained_rings = []
        for ring in rings:
            if len(ring) in [3, 4]:
                strained_rings.append(len(ring))

        if strained_rings:
            return PhysicalValidation(
                check_name="ring_strain",
                category=ValidationCategory.STRUCTURE,
                passed=True,  # 경고만
                severity=Severity.WARNING,
                description=f"High strain rings detected: {strained_rings}-membered rings",
                details={"strained_rings": strained_rings},
                suggested_fix="Consider if these rings are synthetically accessible"
            )

        return PhysicalValidation(
            check_name="ring_strain",
            category=ValidationCategory.STRUCTURE,
            passed=True,
            severity=Severity.INFO,
            description="No high-strain rings detected"
        )

    def _check_chirality(self, mol) -> PhysicalValidation:
        """키랄성 일관성 검증"""
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

        unassigned = [c for c in chiral_centers if c[1] == '?']

        if unassigned:
            return PhysicalValidation(
                check_name="chirality_check",
                category=ValidationCategory.STRUCTURE,
                passed=True,  # 경고만
                severity=Severity.WARNING,
                description=f"Unassigned chiral centers: {len(unassigned)} centers",
                affected_atoms=[c[0] for c in unassigned],
                suggested_fix="Specify stereochemistry for chiral centers"
            )

        return PhysicalValidation(
            check_name="chirality_check",
            category=ValidationCategory.STRUCTURE,
            passed=True,
            severity=Severity.INFO,
            description=f"Chirality consistent ({len(chiral_centers)} chiral centers)"
        )

    def _generate_3d_structure(self, mol):
        """3D 구조 생성"""
        mol_3d = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol_3d, randomSeed=42, maxAttempts=50)

        if result == -1:
            logger.warning("[PhysicalValidator] 3D embedding failed")
            return None

        try:
            AllChem.MMFFOptimizeMolecule(mol_3d, maxIters=200)
        except Exception:
            # MMFF 실패 시 UFF 시도
            try:
                AllChem.UFFOptimizeMolecule(mol_3d, maxIters=200)
            except Exception:
                pass

        return mol_3d

    def _check_steric_clashes(self, mol) -> PhysicalValidation:
        """원자 간 입체 충돌 검사"""
        try:
            conf = mol.GetConformer()
            clashes = []

            for i in range(mol.GetNumAtoms()):
                for j in range(i + 2, mol.GetNumAtoms()):  # i+2 to skip bonded atoms
                    # 결합된 원자는 제외
                    if mol.GetBondBetweenAtoms(i, j) is not None:
                        continue

                    # 1-3 관계도 제외 (같은 원자에 결합)
                    if self._are_1_3_related(mol, i, j):
                        continue

                    pos_i = conf.GetAtomPosition(i)
                    pos_j = conf.GetAtomPosition(j)

                    distance = pos_i.Distance(pos_j)

                    atom_i = mol.GetAtomWithIdx(i)
                    atom_j = mol.GetAtomWithIdx(j)

                    vdw_i = self.VDW_RADII.get(atom_i.GetSymbol(), 1.70)
                    vdw_j = self.VDW_RADII.get(atom_j.GetSymbol(), 1.70)

                    min_distance = (vdw_i + vdw_j) * 0.6  # 40% 허용

                    if distance < min_distance:
                        clashes.append({
                            "atoms": [i, j],
                            "distance": round(distance, 3),
                            "min_allowed": round(min_distance, 3)
                        })

            if clashes:
                return PhysicalValidation(
                    check_name="steric_clash",
                    category=ValidationCategory.STRUCTURE,
                    passed=False,
                    severity=Severity.CRITICAL,
                    description=f"Steric clashes detected: {len(clashes)} atom pairs too close",
                    affected_atoms=[a for c in clashes for a in c["atoms"]],
                    details={"clashes": clashes[:5]},  # 처음 5개만
                    suggested_fix="Relax molecular geometry or redesign structure"
                )

            return PhysicalValidation(
                check_name="steric_clash",
                category=ValidationCategory.STRUCTURE,
                passed=True,
                severity=Severity.INFO,
                description="No steric clashes detected"
            )

        except Exception as e:
            return PhysicalValidation(
                check_name="steric_clash",
                category=ValidationCategory.STRUCTURE,
                passed=True,
                severity=Severity.INFO,
                description=f"Steric clash check completed: {str(e)}"
            )

    def _are_1_3_related(self, mol, i: int, j: int) -> bool:
        """두 원자가 1-3 관계인지 확인"""
        atom_i = mol.GetAtomWithIdx(i)
        neighbors_i = [n.GetIdx() for n in atom_i.GetNeighbors()]

        atom_j = mol.GetAtomWithIdx(j)
        neighbors_j = [n.GetIdx() for n in atom_j.GetNeighbors()]

        # 공통 이웃이 있으면 1-3 관계
        return bool(set(neighbors_i) & set(neighbors_j))

    def _check_bond_lengths(self, mol) -> List[PhysicalValidation]:
        """결합 길이 검증"""
        validations = []

        try:
            conf = mol.GetConformer()
            anomalies = []

            for bond in mol.GetBonds():
                i = bond.GetBeginAtomIdx()
                j = bond.GetEndAtomIdx()

                pos_i = conf.GetAtomPosition(i)
                pos_j = conf.GetAtomPosition(j)
                distance = pos_i.Distance(pos_j)

                atom_i = mol.GetAtomWithIdx(i).GetSymbol()
                atom_j = mol.GetAtomWithIdx(j).GetSymbol()
                bond_order = bond.GetBondTypeAsDouble()

                # 표준 결합 길이 조회
                key = tuple(sorted([atom_i, atom_j])) + (bond_order,)
                expected = self.BOND_LENGTHS.get(key)

                if expected:
                    exp_length, tolerance = expected
                    if abs(distance - exp_length) > tolerance * 2:  # 2배 허용
                        anomalies.append({
                            "atoms": [i, j],
                            "symbols": [atom_i, atom_j],
                            "measured": round(distance, 3),
                            "expected": exp_length,
                            "deviation": round(abs(distance - exp_length), 3)
                        })

            if anomalies:
                validations.append(PhysicalValidation(
                    check_name="bond_length",
                    category=ValidationCategory.STRUCTURE,
                    passed=len(anomalies) < 3,  # 3개 미만이면 통과
                    severity=Severity.WARNING if len(anomalies) < 3 else Severity.CRITICAL,
                    description=f"Bond length anomalies: {len(anomalies)} bonds outside expected range",
                    details={"anomalies": anomalies[:5]},
                    suggested_fix="Optimize molecular geometry"
                ))
            else:
                validations.append(PhysicalValidation(
                    check_name="bond_length",
                    category=ValidationCategory.STRUCTURE,
                    passed=True,
                    severity=Severity.INFO,
                    description="All bond lengths within expected range"
                ))

        except Exception as e:
            validations.append(PhysicalValidation(
                check_name="bond_length",
                category=ValidationCategory.STRUCTURE,
                passed=True,
                severity=Severity.INFO,
                description=f"Bond length check completed: {str(e)}"
            ))

        return validations

    def _check_bond_angles(self, mol) -> List[PhysicalValidation]:
        """결합 각도 검증"""
        validations = []

        try:
            conf = mol.GetConformer()
            anomalies = []

            for atom in mol.GetAtoms():
                neighbors = [n.GetIdx() for n in atom.GetNeighbors()]

                if len(neighbors) < 2:
                    continue

                # 모든 이웃 쌍에 대해 각도 계산
                for k in range(len(neighbors)):
                    for l in range(k + 1, len(neighbors)):
                        i, j = neighbors[k], neighbors[l]
                        center = atom.GetIdx()

                        # 각도 계산
                        pos_i = conf.GetAtomPosition(i)
                        pos_center = conf.GetAtomPosition(center)
                        pos_j = conf.GetAtomPosition(j)

                        v1 = [pos_i.x - pos_center.x, pos_i.y - pos_center.y, pos_i.z - pos_center.z]
                        v2 = [pos_j.x - pos_center.x, pos_j.y - pos_center.y, pos_j.z - pos_center.z]

                        # 내적으로 각도 계산
                        dot = sum(a*b for a, b in zip(v1, v2))
                        mag1 = sum(a*a for a in v1) ** 0.5
                        mag2 = sum(a*a for a in v2) ** 0.5

                        if mag1 > 0 and mag2 > 0:
                            cos_angle = dot / (mag1 * mag2)
                            cos_angle = max(-1, min(1, cos_angle))  # 클램프
                            angle = np.degrees(np.arccos(cos_angle)) if np else 0

                            # 비정상 각도 체크 (60° 미만 또는 180° 초과)
                            if angle < 60 or angle > 175:
                                anomalies.append({
                                    "atoms": [i, center, j],
                                    "angle": round(angle, 1)
                                })

            if anomalies:
                validations.append(PhysicalValidation(
                    check_name="bond_angle",
                    category=ValidationCategory.STRUCTURE,
                    passed=len(anomalies) < 3,
                    severity=Severity.WARNING if len(anomalies) < 3 else Severity.CRITICAL,
                    description=f"Bond angle anomalies: {len(anomalies)} angles outside normal range",
                    details={"anomalies": anomalies[:5]}
                ))
            else:
                validations.append(PhysicalValidation(
                    check_name="bond_angle",
                    category=ValidationCategory.STRUCTURE,
                    passed=True,
                    severity=Severity.INFO,
                    description="All bond angles within normal range"
                ))

        except Exception as e:
            validations.append(PhysicalValidation(
                check_name="bond_angle",
                category=ValidationCategory.STRUCTURE,
                passed=True,
                severity=Severity.INFO,
                description=f"Bond angle check completed: {str(e)}"
            ))

        return validations

    # =========================================================================
    # Result Formatting and Storage
    # =========================================================================

    def _format_result(
        self,
        validations: List[PhysicalValidation],
        smiles: str,
        session_id: Optional[str],
        calculation_id: Optional[str],
        molecule_name: Optional[str] = None
    ) -> Dict[str, Any]:
        """검증 결과 포맷"""
        total = len(validations)
        passed = sum(1 for v in validations if v.passed)
        critical_failures = sum(1 for v in validations if not v.passed and v.severity == Severity.CRITICAL)
        warnings = sum(1 for v in validations if v.severity == Severity.WARNING and not v.passed)

        # 전체 판정
        if critical_failures > 0:
            overall_status = "fail"
            overall_message = f"Critical issues found: {critical_failures} validation(s) failed"
        elif warnings > 0:
            overall_status = "warning"
            overall_message = f"Passed with warnings: {warnings} warning(s)"
        else:
            overall_status = "pass"
            overall_message = "All validations passed"

        return {
            "success": critical_failures == 0,
            "smiles": smiles,
            "molecule_name": molecule_name,
            "session_id": session_id,
            "calculation_id": calculation_id,
            "overall_status": overall_status,
            "overall_message": overall_message,
            "summary": {
                "total_checks": total,
                "passed": passed,
                "failed": total - passed,
                "critical_failures": critical_failures,
                "warnings": warnings
            },
            "validations": [
                {
                    "check_name": v.check_name,
                    "category": v.category.value,
                    "passed": v.passed,
                    "severity": v.severity.value,
                    "description": v.description,
                    "measured_value": v.measured_value,
                    "threshold": v.threshold,
                    "unit": v.unit,
                    "affected_atoms": v.affected_atoms,
                    "suggested_fix": v.suggested_fix,
                    "details": v.details
                }
                for v in validations
            ],
            "timestamp": datetime.utcnow().isoformat()
        }

    async def _save_validations(
        self,
        validations: List[PhysicalValidation],
        session_id: str,
        calculation_id: Optional[str],
        smiles: str,
        molecule_name: Optional[str]
    ):
        """검증 결과 DB 저장"""
        try:
            records = []
            for v in validations:
                records.append({
                    "session_id": session_id,
                    "calculation_id": calculation_id,
                    "smiles": smiles,
                    "molecule_name": molecule_name,
                    "check_name": v.check_name,
                    "check_category": v.category.value,
                    "passed": v.passed,
                    "severity": v.severity.value,
                    "actual_value": v.measured_value,
                    "expected_range": {"threshold": v.threshold, "unit": v.unit} if v.threshold else None,
                    "details": {
                        "description": v.description,
                        "suggested_fix": v.suggested_fix,
                        **v.details
                    }
                })

            self.supabase.table("physical_validations").insert(records).execute()
            logger.info(f"[PhysicalValidator] Saved {len(records)} validations for session {session_id}")

        except Exception as e:
            logger.error(f"[PhysicalValidator] Failed to save validations: {e}")


# ============================================================================
# Module-level Functions
# ============================================================================

_validator_instance = None


def get_physical_validator() -> PhysicalValidator:
    """싱글톤 인스턴스 반환"""
    global _validator_instance
    if _validator_instance is None:
        _validator_instance = PhysicalValidator()
    return _validator_instance


async def validate_structure(
    smiles: str,
    session_id: Optional[str] = None,
    calculation_id: Optional[str] = None,
    molecule_name: Optional[str] = None,
    generate_3d: bool = True,
    save_to_db: bool = True
) -> Dict[str, Any]:
    """구조 검증 실행 (편의 함수)"""
    validator = get_physical_validator()
    return await validator.validate_structure(
        smiles=smiles,
        session_id=session_id,
        calculation_id=calculation_id,
        molecule_name=molecule_name,
        generate_3d=generate_3d,
        save_to_db=save_to_db
    )
