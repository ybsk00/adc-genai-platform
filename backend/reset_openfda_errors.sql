-- OpenFDA 오류 데이터 초기화
-- relevance_score = 0 또는 target이 null인 OpenFDA 데이터를 재분석 대기 상태로 되돌림

UPDATE golden_set_library 
SET 
    ai_refined = false, 
    relevance_score = 0, 
    outcome_type = 'Unknown',  -- 'Pending'은 체크 제약 위반, 'Unknown' 사용
    properties = properties - 'ai_analysis'
WHERE enrichment_source = 'open_fda_api' 
  AND (
    properties->'ai_analysis'->>'target' IS NULL 
    OR relevance_score = 0
    OR ai_refined = true
  );

-- 결과 확인
SELECT 
    COUNT(*) as reset_count,
    COUNT(*) FILTER (WHERE ai_refined = false) as pending_count
FROM golden_set_library 
WHERE enrichment_source = 'open_fda_api';
