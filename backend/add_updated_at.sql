-- Add updated_at column if it doesn't exist
ALTER TABLE golden_set_library 
ADD COLUMN IF NOT EXISTS updated_at TIMESTAMPTZ DEFAULT NOW();

-- Create extension for automatic timestamp updates if not exists
CREATE EXTENSION IF NOT EXISTS moddatetime;

-- Create trigger to automatically update updated_at
DROP TRIGGER IF EXISTS handle_updated_at ON golden_set_library;

CREATE TRIGGER handle_updated_at
BEFORE UPDATE ON golden_set_library
FOR EACH ROW
EXECUTE PROCEDURE moddatetime (updated_at);
