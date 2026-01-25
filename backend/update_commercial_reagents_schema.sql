-- commercial_reagents 테이블에 정량 지표 및 수동 수정 플래그 추가
ALTER TABLE commercial_reagents 
ADD COLUMN IF NOT EXISTS is_manual_override BOOLEAN DEFAULT FALSE,
ADD COLUMN IF NOT EXISTS binding_affinity TEXT,
ADD COLUMN IF NOT EXISTS isotype TEXT,
ADD COLUMN IF NOT EXISTS host_species TEXT,
ADD COLUMN IF NOT EXISTS orr_pct TEXT,
ADD COLUMN IF NOT EXISTS os_months TEXT,
ADD COLUMN IF NOT EXISTS pfs_months TEXT;

-- 기존 데이터의 ai_refined가 null인 경우 false로 초기화
UPDATE commercial_reagents SET ai_refined = FALSE WHERE ai_refined IS NULL;