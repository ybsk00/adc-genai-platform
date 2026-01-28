# ADC 설계 플랫폼 UI 고도화 구현 계획

**작성일**: 2026-01-28
**기준 문서**: UI_UX개선.md (v2.3)
**현재 상태**: Phase 3 구현 완료 (UniProt API + Environment Lock)

---

## 1. 구현 완료 항목 ✅

### 1.1 공통 레이아웃 (Trinity Layout)
- ✅ `TrinityLayout.tsx` - 데스크톱용 3단 분할 레이아웃
- ✅ `TrinityTabsLayout` - 모바일용 탭 방식 레이아웃

### 1.2 시각화 컴포넌트
- ✅ `MoleculeViewer2D.tsx` - smiles-drawer 기반 2D 구조 뷰어
- ✅ `MoleculeDiffViewer` - 두 분자 나란히 비교
- ✅ `SuccessRadarChart.tsx` - Recharts 기반 거미줄 차트
- ✅ `SAScoreGauge.tsx` - 합성 용이성 1~10점 게이지

### 1.3 Evidence Hub
- ✅ `EvidenceHub.tsx` - Librarian 근거 논문 패널
- ✅ 문헌 참조 리스트 (PMID 연동)
- ✅ Golden Set 참조 리스트
- ✅ 핀포인트 근거 하이라이트 (스캐폴드-논문 연동)

### 1.4 Audit & Compliance
- ✅ `DigitalSealBadge.tsx` - SHA-256 해시 배지
- ✅ `CompactSealBadge` - 인라인 사용 배지
- ✅ `AuditorFeedback.tsx` - Redesign Loop 상태 시각화
- ✅ Audit Checklist 컴포넌트 (Pass/Fail)
- ✅ Reasoning Popover (추론 근거)

### 1.5 Optimization Metrics
- ✅ `DeltaMetricDisplay.tsx` - 전후 비교 화살표 지표
- ✅ `CompactDelta` - 인라인 델타 표시

### 1.6 Self-Healing Visualization
- ✅ `HealerAnimation.tsx` - 자가 치유 라이브 중계
- ✅ 코드 에러 수정 애니메이션 (타이핑 효과)
- ✅ `HealerStatusIndicator` - 상태 표시 배지

### 1.7 CMC & Sourcing
- ✅ `BOMTable.tsx` - Bill of Materials 테이블
- ✅ `BOMSummary` - 요약 뱃지

### 1.8 WebSocket 이벤트 (백엔드)
- ✅ `shared_state_sync` - SMILES/메트릭 동기화
- ✅ `healer_action` - 자가 치유 상태
- ✅ `librarian_references` - 논문/Golden Set 참조
- ✅ `auditor_feedback` - Redesign Loop 피드백
- ✅ `digital_seal` - 해시 정보

### 1.9 백엔드 API
- ✅ `/api/design/session/{id}/digital-seal` - 해시 조회
- ✅ `/api/design/session/{id}/audit-bundle` - Audit Trail 내보내기

### 1.10 Hook 업데이트
- ✅ `useDesignSession.ts` - 모든 새 이벤트 핸들러 추가

### 1.11 페이지 업데이트 (Phase 2 완료)
- ✅ `DenovoDesign.tsx` - Trinity Layout 적용
- ✅ `LeadOptimization.tsx` - Trinity + MoleculeDiffViewer + DeltaMetricDisplay
- ✅ `PreclinicalAudit.tsx` - Trinity + AuditorFeedback + DigitalSealBadge
- ✅ `CMCSourcing.tsx` - Trinity + BOMTable + SAScoreGauge

---

## 2. 구현된 사용자 요청 기능

### 2.1 Librarian '핀포인트 근거' 시각화
- 중앙 Molecule Viewer에서 스캐폴드 클릭 → 우측 Evidence Hub에서 관련 논문 하이라이트
- `selectedScaffoldId` 상태로 연동

### 2.2 Auditor '반려 및 피드백' UI
- `AuditorFeedback.tsx`의 `RedesignLoopAlert` 컴포넌트
- "DAR 목표 미달로 인해 Alchemist에게 구조 변경 요청 중" 메시지 노출
- 반복 횟수 (Iteration N/3) 표시

### 2.3 'Audit Trail Bundle' 내보내기
- `DigitalSealBadge`의 다운로드 버튼
- `/api/design/session/{id}/audit-bundle` API
- 포함 항목: SMILES, 코드, 추론 로그, SHA-256 해시

### 2.4 Healer '자가 치유 라이브 중계'
- `HealerAnimation.tsx` - 애니메이션 효과
- 빨간색 에러 밑줄 → 코드 타이핑 애니메이션
- 상태별 배지 (Detecting → Analyzing → Healing → Healed)

---

## 3. 파일 구조

```
frontend/src/
├── components/design/
│   ├── index.ts                    # ✅ Export index
│   ├── TrinityLayout.tsx           # ✅ 3단 분할 레이아웃
│   ├── MoleculeViewer2D.tsx        # ✅ 2D 분자 뷰어
│   ├── SuccessRadarChart.tsx       # ✅ Radar 차트
│   ├── EvidenceHub.tsx             # ✅ Librarian 근거
│   ├── DigitalSealBadge.tsx        # ✅ 21 CFR Part 11 배지
│   ├── DeltaMetricDisplay.tsx      # ✅ 전후 비교
│   ├── HealerAnimation.tsx         # ✅ 자가 치유 애니메이션
│   ├── AuditorFeedback.tsx         # ✅ Redesign Loop
│   ├── BOMTable.tsx                # ✅ Bill of Materials
│   ├── SAScoreGauge.tsx            # ✅ SA Score 게이지
│   ├── AntigenAutocomplete.tsx     # ✅ UniProt 자동완성 (Phase 3)
│   ├── EnvironmentLock.tsx         # ✅ 환경 버전 표시 (Phase 3)
│   ├── LiveAgentConsole.tsx        # 기존
│   └── DesignProgressTimeline.tsx  # 기존
├── hooks/
│   └── useDesignSession.ts         # ✅ 업데이트됨
└── pages/Dashboard/
    ├── DenovoDesign.tsx            # ✅ Trinity Layout 적용
    ├── LeadOptimization.tsx        # ✅ Trinity Layout 적용 (Phase 2)
    ├── PreclinicalAudit.tsx        # ✅ Trinity Layout 적용 (Phase 2)
    └── CMCSourcing.tsx             # ✅ Trinity Layout 적용 (Phase 2)

backend/app/
├── core/
│   └── websocket_hub.py            # ✅ 새 이벤트 추가
└── api/
    ├── design.py                   # ✅ 새 API 추가
    ├── uniprot.py                  # ✅ UniProt 검색 API (Phase 3)
    └── system.py                   # ✅ 환경 정보 API (Phase 3)
```

---

## 4. 결정 사항 (확정)

| 항목 | 결정 |
|------|------|
| 분자 구조 뷰어 | **smiles-drawer** (RDKit.js WASM은 추후 업그레이드) |
| Trinity Layout | **4대 메뉴 전체 적용** |
| 백엔드 API | **Digital Seal 해시 조회 API 1순위** |

---

## 5. Phase 2 완료 ✅

### 5.1 모든 페이지 Trinity Layout 적용 완료
- ✅ LeadOptimization.tsx - Trinity + MoleculeDiffViewer + DeltaMetricDisplay
- ✅ PreclinicalAudit.tsx - Trinity + AuditorFeedback + DigitalSealBadge
- ✅ CMCSourcing.tsx - Trinity + BOMTable + SAScoreGauge

---

## 6. Phase 3 구현 완료 ✅

### 6.1 추가 기능 (완료)
- ✅ UniProt 항원 검색 API
  - `backend/app/api/uniprot.py` - 검색, 상세조회, 자주 사용 타겟 목록
  - `frontend/src/components/design/AntigenAutocomplete.tsx` - 자동완성 컴포넌트
  - DenovoDesign, PreclinicalAudit 페이지에 적용
- ✅ Environment Lock 버전 표시
  - `backend/app/api/system.py` - 환경 정보, 버전 잠금 API
  - `frontend/src/components/design/EnvironmentLock.tsx` - 버전 표시 컴포넌트
  - 모든 4대 메뉴 Right Panel에 적용

### 6.2 Premium 기능 (P3)
- [ ] 3D Viewer (Mol*) - AlphaFold 연동
- [x] RDKit.js WASM 업그레이드 ✅ (2026-01-28 완료)
  - MoleculeViewer2D.tsx: smiles-drawer → @rdkit/rdkit WASM
  - 새 기능: SMARTS 기반 서브구조 하이라이팅, 분자 검증
  - useRDKit 훅 추가: validateSmiles, canonicalizeSmiles, getMolBlock, hasSubstructMatch

---

## 7. 사용법

### 컴포넌트 Import
```tsx
import {
  TrinityTabsLayout,
  MoleculeViewer2D,
  SuccessRadarChart,
  EvidenceHub,
  DigitalSealBadge,
  HealerAnimation,
  AuditorFeedback,
  DeltaMetricDisplay,
  BOMTable,
  SAScoreGauge
} from '@/components/design'
```

### Hook 사용
```tsx
const {
  isConnected,
  agentLogs,
  currentSmiles,
  calculatedMetrics,
  scaffolds,
  healerAction,
  isHealerActive,
  literatureRefs,
  goldenSetRefs,
  auditFeedback,
  digitalSeal
} = useDesignSession({ sessionId })
```

---

*Phase 1, 2, 3 구현 완료.*
*- Phase 1: 모든 UI 컴포넌트 및 WebSocket 이벤트*
*- Phase 2: 4대 메뉴 Trinity Layout 적용*
*- Phase 3: UniProt API 자동완성, Environment Lock 버전 표시*
*다음 단계: 3D Viewer (Mol*), RDKit.js WASM 등 Premium 기능*
