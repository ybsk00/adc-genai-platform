-- Add mdl_number column to commercial_reagents table
ALTER TABLE commercial_reagents ADD COLUMN IF NOT EXISTS mdl_number TEXT;
