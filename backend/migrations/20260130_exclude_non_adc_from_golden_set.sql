-- ============================================================================
-- Golden Set 데이터 정제: 비-ADC 약물 제외 및 이름 정규화
-- 실행일: 2026-01-30
-- ============================================================================

-- ── 1. 제외 대상 확인 (DRY RUN) ──
-- 먼저 어떤 데이터가 제외되는지 확인하세요
SELECT id, name, target_1, category, status
FROM golden_set_library
WHERE name ILIKE ANY(ARRAY[
    '%Herceptin%',
    '%Trazimera%',
    '%OGIVRI%',
    '%Kanjinti%',
    '%Phesgo%',
    '%Herceptin Hylecta%'
]);

-- ── 2. 비-ADC 약물 상태 변경 (soft delete) ──
-- status를 'rejected'로 변경하고 reviewer_note에 사유 기록
UPDATE golden_set_library
SET
    status = 'rejected',
    reviewer_note = '비-ADC 약물 제외: 순수 단일클론항체(mAb) 또는 항체 복합제. antibody_library로 이동 대상.',
    review_required = false,
    updated_at = now()
WHERE name ILIKE ANY(ARRAY[
    '%Herceptin%',
    '%Trazimera%',
    '%OGIVRI%',
    '%Kanjinti%',
    '%Phesgo%'
]);

-- ── 3. 오타 교정 및 이름 정규화 ──
-- DATROWAY → Dato-DXd (Datopotamab deruxtecan)
UPDATE golden_set_library
SET
    name = 'Dato-DXd (Datopotamab deruxtecan)',
    description = COALESCE(description, '') || ' [정규화: DATROWAY → Dato-DXd]',
    updated_at = now()
WHERE name ILIKE '%DATROWAY%';

-- EMRELIS → Telisotuzumab vedotin
UPDATE golden_set_library
SET
    name = 'Emrelis (Telisotuzumab vedotin)',
    description = COALESCE(description, '') || ' [정규화: EMRELIS → Telisotuzumab vedotin]',
    updated_at = now()
WHERE name ILIKE '%EMRELIS%';

-- ── 4. 정제 결과 확인 ──
-- 제외된 항목
SELECT name, status, reviewer_note
FROM golden_set_library
WHERE status = 'rejected'
  AND reviewer_note ILIKE '%비-ADC%'
ORDER BY name;

-- 유지되는 ADC 목록
SELECT name, target_1, target_2, status, orr_pct, pfs_months, os_months
FROM golden_set_library
WHERE status != 'rejected'
  AND category ILIKE '%ADC%' OR name ILIKE ANY(ARRAY[
    '%KADCYLA%', '%TRODELVY%', '%Enhertu%', '%Dato-DXd%',
    '%TIVDAK%', '%Emrelis%', '%PADCEV%', '%ADCETRIS%', '%POLIVY%'
])
ORDER BY orr_pct DESC NULLS LAST;

-- ── 5. (선택) 완전 삭제가 필요한 경우 ──
-- 주의: 이 쿼리는 데이터를 영구 삭제합니다
-- DELETE FROM golden_set_library
-- WHERE status = 'rejected'
--   AND reviewer_note ILIKE '%비-ADC%';
