-- Add missing columns to commercial_reagents table
ALTER TABLE commercial_reagents ADD COLUMN IF NOT EXISTS target text;
ALTER TABLE commercial_reagents ADD COLUMN IF NOT EXISTS properties jsonb;
ALTER TABLE commercial_reagents ADD COLUMN IF NOT EXISTS source_name text;
ALTER TABLE commercial_reagents ADD COLUMN IF NOT EXISTS summary text;
ALTER TABLE commercial_reagents ADD COLUMN IF NOT EXISTS price_data jsonb;

-- Optional: Add index for source_name if needed
CREATE INDEX IF NOT EXISTS idx_commercial_reagents_source_name ON commercial_reagents(source_name);
