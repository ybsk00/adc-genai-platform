-- [2026-01-21] outcome_type 제약 조건 업데이트
-- 기존 제약 조건 삭제 후 'Ongoing'을 포함하여 재생성

ALTER TABLE public.golden_set_library 
DROP CONSTRAINT IF EXISTS golden_set_library_outcome_type_check;

ALTER TABLE public.golden_set_library 
ADD CONSTRAINT golden_set_library_outcome_type_check 
CHECK (outcome_type IN ('Success', 'Failure', 'Terminated', 'Ongoing', 'Unknown'));

-- 확인용
SELECT constraint_name, check_clause 
FROM information_schema.check_constraints 
WHERE constraint_name = 'golden_set_library_outcome_type_check';
