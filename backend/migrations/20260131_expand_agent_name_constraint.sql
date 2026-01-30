-- ============================================================================
-- Migration: Expand agent_name_check constraint
-- Date: 2026-01-31
-- Purpose: Allow all Navigator agent names in agent_execution_logs
-- ============================================================================

-- Drop the existing restrictive constraint
ALTER TABLE agent_execution_logs
DROP CONSTRAINT IF EXISTS agent_name_check;

-- Re-create with expanded agent names (including orchestrator, virtual_trial, etc.)
ALTER TABLE agent_execution_logs
ADD CONSTRAINT agent_name_check CHECK (
    agent_name IN (
        'librarian',
        'alchemist',
        'coder',
        'healer',
        'auditor',
        'orchestrator',
        'virtual_trial',
        'clinical',
        'navigator',
        'competitor'
    )
);

-- Also drop the foreign key constraint that requires session_id to exist in design_sessions
-- (Navigator uses navigator_sessions, not design_sessions)
ALTER TABLE agent_execution_logs
DROP CONSTRAINT IF EXISTS agent_execution_logs_session_id_fkey;

-- Add critical_errors column to navigator_sessions (if not exists)
DO $$
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM information_schema.columns
        WHERE table_name = 'navigator_sessions'
        AND column_name = 'critical_errors'
    ) THEN
        ALTER TABLE navigator_sessions
        ADD COLUMN critical_errors JSONB DEFAULT '[]'::jsonb;
    END IF;
END $$;

-- Add comment for documentation
COMMENT ON COLUMN navigator_sessions.critical_errors IS
    'Fail-Fast: Critical errors that prevented complete computation. No fake data was generated.';
