-- ============================================================================
-- 1. HERZUMA, Ontruzant 제거 (허셉틴 바이오시밀러 - ADC 아님)
-- ============================================================================

-- 확인
SELECT id, name, target_1, status, category
FROM golden_set_library
WHERE name ILIKE ANY(ARRAY['%HERZUMA%', '%Ontruzant%']);

-- 삭제 (soft delete)
UPDATE golden_set_library
SET
    status = 'rejected',
    reviewer_note = '바이오시밀러(Trastuzumab biosimilar) - ADC 아님. 골든셋에서 제외.',
    review_required = false,
    updated_at = now()
WHERE name ILIKE ANY(ARRAY['%HERZUMA%', '%Ontruzant%']);

-- knowledge_base에도 있으면 확인
SELECT id, title, source_type, mapped_golden_set_id
FROM knowledge_base
WHERE title ILIKE ANY(ARRAY['%HERZUMA%', '%Ontruzant%']);


-- ============================================================================
-- 2. 실패 ADC 사례 등록: Vadastuximab talirine (SGN-CD33A)
-- ============================================================================

-- 이미 있는지 확인
SELECT id, name, target_1, outcome_type, status
FROM golden_set_library
WHERE name ILIKE ANY(ARRAY['%SGN-CD33A%', '%vadastuximab%', '%바다스투맙%']);

-- 없으면 INSERT, 있으면 UPDATE
INSERT INTO golden_set_library (
    name, category, description, target_1, target_symbol,
    outcome_type, failure_reason, status,
    properties, orr_pct, pfs_months, os_months
) VALUES (
    'Vadastuximab talirine (SGN-CD33A)',
    'ADC / Hematologic Malignancy / AML',
    'CD33 타겟 ADC. 급성골수성백혈병(AML) 치료 목적. '
    '임상 중 치명적 간 정맥 폐쇄 질환(VOD/SOS) 발생으로 개발 중단. '
    'AI 학습 포인트: 간독성 유발 구조 회피, Safety 모델 강화.',
    'CD33',
    'CD33',
    'Failure',
    '간독성(hepatotoxicity): 임상 중 간 정맥 폐쇄 질환(VOD/SOS) 등 치명적 간독성 발생으로 임상 중단.',
    'approved',
    jsonb_build_object(
        'payload_class', 'PBD dimer (SGD-1882)',
        'linker_type', 'cleavable',
        'mechanism_of_action', 'DNA crosslinker',
        'dar', 2,
        'learning_point', 'Safety 모델 강화: PBD dimer 페이로드의 간독성 위험을 AI가 경고해야 함.',
        'failure_category', 'safety_toxicity',
        'clinical_phase_at_failure', 'Phase 3'
    ),
    NULL, NULL, NULL
)
ON CONFLICT DO NOTHING;


-- ============================================================================
-- 3. 실패 ADC 사례 등록: Glembatumumab vedotin (CDX-011)
-- ============================================================================

-- 이미 있는지 확인
SELECT id, name, target_1, outcome_type, status
FROM golden_set_library
WHERE name ILIKE ANY(ARRAY['%CDX-011%', '%glembatumumab%', '%글렘바투무맙%']);

-- 없으면 INSERT, 있으면 UPDATE
INSERT INTO golden_set_library (
    name, category, description, target_1, target_symbol,
    outcome_type, failure_reason, status,
    properties, orr_pct, pfs_months, os_months
) VALUES (
    'Glembatumumab vedotin (CDX-011)',
    'ADC / Breast Cancer / TNBC',
    'gpNMB 타겟 ADC. 삼중음성유방암(TNBC) 대상 METRIC 임상시험에서 '
    '대조군(카페시타빈) 대비 생존 기간 개선 실패. '
    'AI 학습 포인트: 타겟 발현율과 실제 약효 간 상관관계를 냉정하게 판단.',
    'gpNMB',
    'GPNMB',
    'Failure',
    '효능 미달(lack of efficacy): METRIC trial에서 TNBC 대상으로 대조군 대비 OS/PFS 개선 실패.',
    'approved',
    jsonb_build_object(
        'payload_class', 'MMAE',
        'linker_type', 'cleavable (vc)',
        'mechanism_of_action', 'Microtubule inhibitor',
        'dar', 4,
        'learning_point', 'Efficacy 예측 정밀도: 타겟 발현율이 높아도 실제 약효와 상관관계가 낮을 수 있음. AI가 더 냉정하게 판단해야 함.',
        'failure_category', 'efficacy_failure',
        'clinical_phase_at_failure', 'Phase 2',
        'trial_name', 'METRIC trial'
    ),
    NULL, NULL, NULL
)
ON CONFLICT DO NOTHING;


-- ============================================================================
-- 4. 기존에 이미 있었다면 → Failure로 업데이트
-- ============================================================================

-- SGN-CD33A가 이미 있으면 Failure로 변경
UPDATE golden_set_library
SET
    outcome_type = 'Failure',
    failure_reason = '간독성(hepatotoxicity): 임상 중 간 정맥 폐쇄 질환(VOD/SOS) 등 치명적 간독성 발생으로 임상 중단.',
    properties = COALESCE(properties, '{}'::jsonb) || jsonb_build_object(
        'payload_class', 'PBD dimer (SGD-1882)',
        'linker_type', 'cleavable',
        'mechanism_of_action', 'DNA crosslinker',
        'dar', 2,
        'learning_point', 'Safety 모델 강화: PBD dimer 페이로드의 간독성 위험을 AI가 경고해야 함.',
        'failure_category', 'safety_toxicity'
    ),
    updated_at = now()
WHERE name ILIKE ANY(ARRAY['%SGN-CD33A%', '%vadastuximab%'])
  AND outcome_type IS DISTINCT FROM 'Failure';

-- CDX-011이 이미 있으면 Failure로 변경
UPDATE golden_set_library
SET
    outcome_type = 'Failure',
    failure_reason = '효능 미달(lack of efficacy): METRIC trial에서 TNBC 대상으로 대조군 대비 OS/PFS 개선 실패.',
    properties = COALESCE(properties, '{}'::jsonb) || jsonb_build_object(
        'payload_class', 'MMAE',
        'linker_type', 'cleavable (vc)',
        'mechanism_of_action', 'Microtubule inhibitor',
        'dar', 4,
        'learning_point', 'Efficacy 예측 정밀도: 타겟 발현율이 높아도 실제 약효와 상관관계가 낮을 수 있음.',
        'failure_category', 'efficacy_failure',
        'trial_name', 'METRIC trial'
    ),
    updated_at = now()
WHERE name ILIKE ANY(ARRAY['%CDX-011%', '%glembatumumab%'])
  AND outcome_type IS DISTINCT FROM 'Failure';


-- ============================================================================
-- 5. knowledge_base에 실패 사례 학습 데이터 등록
-- ============================================================================

-- SGN-CD33A 실패 학습
INSERT INTO knowledge_base (
    source_type, title, summary, content, relevance_score, source_tier, ai_reasoning
) VALUES (
    'failure_case',
    'ADC 실패 사례: Vadastuximab talirine (SGN-CD33A) - 치명적 간독성',
    'CD33 타겟 ADC. PBD dimer 페이로드 사용 시 간 정맥 폐쇄 질환(VOD/SOS) 발생으로 임상 중단.',
    E'## Vadastuximab talirine (SGN-CD33A) 실패 분석\n\n'
    '### 기본 정보\n'
    '- **타겟**: CD33 (AML 세포 표면 항원)\n'
    '- **페이로드**: PBD dimer (SGD-1882)\n'
    '- **링커**: Cleavable\n'
    '- **DAR**: 2\n\n'
    '### 실패 원인\n'
    '임상 Phase 3 중 치명적 간 정맥 폐쇄 질환(VOD/SOS) 발생.\n'
    'PBD dimer 계열 페이로드의 DNA 교차결합 메커니즘이 간세포에도 비특이적 독성을 유발.\n\n'
    '### AI 학습 포인트\n'
    '1. **PBD dimer 페이로드 사용 시 간독성 위험 경고** 필수\n'
    '2. DAR이 낮아도(2) 고독성 페이로드는 안전성 문제 발생 가능\n'
    '3. CD33 타겟 ADC 설계 시 페이로드 선택에 각별한 주의\n'
    '4. Safety 모델에서 간독성 유발 구조 패턴 학습 필요',
    0.95,
    1,
    'ADC 실패 사례 분석: PBD dimer 페이로드의 간독성 위험을 Safety 모델에 반영해야 함.'
)
ON CONFLICT DO NOTHING;

-- CDX-011 실패 학습
INSERT INTO knowledge_base (
    source_type, title, summary, content, relevance_score, source_tier, ai_reasoning
) VALUES (
    'failure_case',
    'ADC 실패 사례: Glembatumumab vedotin (CDX-011) - TNBC 효능 미달',
    'gpNMB 타겟 ADC. METRIC trial에서 삼중음성유방암(TNBC) 대상으로 대조군(카페시타빈) 대비 생존 기간 개선 실패.',
    E'## Glembatumumab vedotin (CDX-011) 실패 분석\n\n'
    '### 기본 정보\n'
    '- **타겟**: gpNMB (Glycoprotein NMB)\n'
    '- **페이로드**: MMAE\n'
    '- **링커**: Cleavable (vc)\n'
    '- **DAR**: 4\n\n'
    '### 실패 원인\n'
    'METRIC 임상시험에서 삼중음성유방암(TNBC) 환자 대상으로\n'
    '대조군(카페시타빈) 대비 전체 생존 기간(OS), 무진행 생존 기간(PFS) 모두 개선 실패.\n'
    'gpNMB 발현이 높은 환자 하위군에서도 유의미한 차이 없음.\n\n'
    '### AI 학습 포인트\n'
    '1. **타겟 발현율 ≠ 약효**: 높은 타겟 발현이 반드시 임상 효과를 보장하지 않음\n'
    '2. TNBC에서 gpNMB 타겟팅은 충분한 치료 효과를 기대하기 어려움\n'
    '3. Efficacy 예측 모델에서 타겟-질환 상관관계를 더 냉정하게 평가 필요\n'
    '4. MMAE+vc 링커 조합은 검증되었으나 타겟 선정이 핵심 실패 요인',
    0.95,
    1,
    'ADC 실패 사례 분석: 타겟 발현율과 실제 약효 간의 상관관계를 AI가 더 냉정하게 판단하도록 학습.'
)
ON CONFLICT DO NOTHING;


-- ============================================================================
-- 6. 최종 확인
-- ============================================================================

-- 제외된 바이오시밀러
SELECT name, status, reviewer_note
FROM golden_set_library
WHERE name ILIKE ANY(ARRAY['%HERZUMA%', '%Ontruzant%', '%Herceptin%', '%Trazimera%', '%OGIVRI%', '%Kanjinti%', '%Phesgo%']);

-- 실패 ADC 사례
SELECT name, target_1, outcome_type, failure_reason
FROM golden_set_library
WHERE outcome_type = 'Failure'
ORDER BY name;

-- knowledge_base 학습 데이터
SELECT id, title, source_type, relevance_score
FROM knowledge_base
WHERE source_type = 'failure_case'
ORDER BY created_at DESC;
