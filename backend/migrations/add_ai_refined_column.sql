-- Add ai_refined column to commercial_reagents table
ALTER TABLE commercial_reagents 
ADD COLUMN IF NOT EXISTS ai_refined BOOLEAN DEFAULT FALSE;

-- Add index for performance
CREATE INDEX IF NOT EXISTS idx_commercial_reagents_ai_refined ON commercial_reagents(ai_refined);
