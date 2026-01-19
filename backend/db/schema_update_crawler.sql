-- Commercial Reagents Table for Ambeed Crawler
CREATE TABLE IF NOT EXISTS commercial_reagents (
    id uuid DEFAULT gen_random_uuid() PRIMARY KEY,
    ambeed_cat_no text UNIQUE NOT NULL,
    product_name text,
    category text, -- 'Payload', 'Linker', 'Conjugate'
    cas_number text,
    smiles_code text,
    molecular_weight text,
    formula text,
    price_data jsonb, -- Store price tiers as JSON
    stock_status text,
    product_url text,
    crawled_at timestamptz DEFAULT now(),
    embedding vector(1536) -- For RAG
);

-- Index for faster search
CREATE INDEX IF NOT EXISTS idx_commercial_reagents_category ON commercial_reagents(category);
CREATE INDEX IF NOT EXISTS idx_commercial_reagents_cas ON commercial_reagents(cas_number);
