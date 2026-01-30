-- ============================================================================
-- Migration: Add agent_logs, warnings columns to navigator_sessions
-- Date: 2026-01-30
-- Reason: V2 orchestrator saves real-time agent logs and warnings
-- ============================================================================

-- 1. agent_logs 컬럼 추가 (실시간 에이전트 실행 로그)
ALTER TABLE navigator_sessions
ADD COLUMN IF NOT EXISTS agent_logs JSONB DEFAULT '[]'::jsonb;

-- 2. warnings 컬럼 추가 (gpNMB, PBD 등 경고 메시지)
ALTER TABLE navigator_sessions
ADD COLUMN IF NOT EXISTS warnings JSONB DEFAULT '[]'::jsonb;

-- 3. status CHECK 제약조건 업데이트 (waiting_for_selection 추가)
ALTER TABLE navigator_sessions
DROP CONSTRAINT IF EXISTS navigator_sessions_status_check;

ALTER TABLE navigator_sessions
ADD CONSTRAINT navigator_sessions_status_check
CHECK (status IN ('pending', 'running', 'completed', 'failed', 'cancelled', 'waiting_for_selection'));

-- 4. 확인
SELECT column_name, data_type, column_default
FROM information_schema.columns
WHERE table_name = 'navigator_sessions'
  AND column_name IN ('agent_logs', 'warnings')
ORDER BY ordinal_position;
