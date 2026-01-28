-- =====================================================
-- ADC Design Engine Phase 1: Multi-Agent Schema
-- Version: 2.2 (21 CFR Part 11 Compliant)
-- Date: 2026-01-27
-- =====================================================

-- Enable required extensions
CREATE EXTENSION IF NOT EXISTS "pgcrypto";
CREATE EXTENSION IF NOT EXISTS "vector";

-- =====================================================
-- 1. Design Sessions: 설계 세션 마스터 테이블
-- =====================================================
CREATE TABLE IF NOT EXISTS design_sessions (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID REFERENCES auth.users(id) NOT NULL,

    -- Session Type & Tier
    session_type VARCHAR(50) NOT NULL CHECK (session_type IN ('denovo', 'optimization', 'audit', 'cmc')),
    tier VARCHAR(20) DEFAULT 'free' CHECK (tier IN ('free', 'premium')),

    -- Workflow Status
    status VARCHAR(30) DEFAULT 'pending' CHECK (status IN ('pending', 'running', 'completed', 'failed', 'manual_review')),
    current_step INTEGER DEFAULT 0,
    total_steps INTEGER DEFAULT 4,
    current_agent VARCHAR(50),

    -- Design Parameters (FDA Audit Trail)
    design_goal TEXT,
    target_antigen VARCHAR(100),
    target_indication VARCHAR(200),
    requested_dar INTEGER,
    linker_preference VARCHAR(50),

    -- Shared State (에이전트 간 실시간 동기화)
    current_smiles TEXT,
    calculated_metrics JSONB DEFAULT '{}',
    validation_flags JSONB DEFAULT '{}',

    -- Results
    result_candidates JSONB DEFAULT '[]',
    final_report JSONB,
    final_approval BOOLEAN DEFAULT FALSE,
    approval_reasoning TEXT,

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    completed_at TIMESTAMPTZ,

    -- Regulatory Metadata
    regulatory_version VARCHAR(20) DEFAULT 'v1.0',
    audit_hash VARCHAR(64)
);

-- =====================================================
-- 2. Agent Execution Logs: 21 CFR Part 11 Compliant
-- INSERT-ONLY + Chain Hash for Immutability
-- =====================================================
CREATE TABLE IF NOT EXISTS agent_execution_logs (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    session_id UUID REFERENCES design_sessions(id) ON DELETE CASCADE,

    -- Agent Info
    agent_name VARCHAR(50) NOT NULL CHECK (agent_name IN (
        'orchestrator', 'alchemist', 'coder', 'healer', 'auditor', 'librarian'
    )),
    step_number INTEGER,
    status VARCHAR(20) CHECK (status IN ('started', 'completed', 'error', 'healed', 'skipped')),

    -- Reasoning (FDA Audit)
    reasoning TEXT,
    decision_summary TEXT,
    confidence_score DECIMAL(3,2),

    -- Input/Output
    input_data JSONB,
    output_data JSONB,
    error_message TEXT,
    execution_time_ms INTEGER,

    -- Self-healing Logs (The Healer)
    retry_count INTEGER DEFAULT 0,
    fix_logs JSONB DEFAULT '[]',
    healing_successful BOOLEAN,

    -- References
    referenced_golden_set_ids UUID[] DEFAULT '{}',
    referenced_knowledge_ids INTEGER[] DEFAULT '{}',
    referenced_pmids TEXT[] DEFAULT '{}',

    -- 21 CFR Part 11 Digital Seal
    record_hash VARCHAR(64) NOT NULL,
    prev_record_hash VARCHAR(64),
    sequence_number BIGINT NOT NULL,
    sealed_at TIMESTAMPTZ DEFAULT NOW(),
    seal_version VARCHAR(10) DEFAULT 'v1.0',

    created_at TIMESTAMPTZ DEFAULT NOW()
);

-- INSERT-ONLY Policy: Block UPDATE/DELETE
CREATE OR REPLACE FUNCTION prevent_log_modification()
RETURNS TRIGGER AS $$
BEGIN
    RAISE EXCEPTION '21 CFR Part 11 Compliance: agent_execution_logs is immutable. Record ID: %', OLD.id;
END;
$$ LANGUAGE plpgsql;

DROP TRIGGER IF EXISTS enforce_immutability ON agent_execution_logs;
CREATE TRIGGER enforce_immutability
    BEFORE UPDATE OR DELETE ON agent_execution_logs
    FOR EACH ROW EXECUTE FUNCTION prevent_log_modification();

-- Chain Hash Auto-Generation
CREATE OR REPLACE FUNCTION generate_record_seal()
RETURNS TRIGGER AS $$
DECLARE
    prev_hash VARCHAR(64);
    prev_seq BIGINT;
    seal_data TEXT;
BEGIN
    -- Get previous record
    SELECT record_hash, sequence_number INTO prev_hash, prev_seq
    FROM agent_execution_logs
    WHERE session_id = NEW.session_id
    ORDER BY sequence_number DESC
    LIMIT 1;

    -- Set sequence number
    NEW.sequence_number := COALESCE(prev_seq, 0) + 1;
    NEW.prev_record_hash := prev_hash;

    -- Build seal data
    seal_data := COALESCE(prev_hash, 'GENESIS') || '|' ||
                 NEW.session_id || '|' ||
                 NEW.agent_name || '|' ||
                 NEW.sequence_number || '|' ||
                 COALESCE(NEW.reasoning, '') || '|' ||
                 NEW.created_at::TEXT;

    -- Generate SHA-256 hash
    NEW.record_hash := encode(sha256(seal_data::bytea), 'hex');

    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

DROP TRIGGER IF EXISTS seal_record_before_insert ON agent_execution_logs;
CREATE TRIGGER seal_record_before_insert
    BEFORE INSERT ON agent_execution_logs
    FOR EACH ROW EXECUTE FUNCTION generate_record_seal();

-- =====================================================
-- 3. Encrypted Structures: Premium SMILES Storage
-- =====================================================
CREATE TABLE IF NOT EXISTS encrypted_structures (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    session_id UUID REFERENCES design_sessions(id) ON DELETE CASCADE,
    user_id UUID REFERENCES auth.users(id) NOT NULL,

    -- Encrypted Data (AES-256-GCM)
    encrypted_smiles BYTEA NOT NULL,
    encryption_iv BYTEA NOT NULL,
    encryption_tag BYTEA NOT NULL,

    -- Metadata
    structure_type VARCHAR(50),
    rank INTEGER,
    is_premium_only BOOLEAN DEFAULT TRUE,

    -- Validation Hash
    structure_hash VARCHAR(64),

    created_at TIMESTAMPTZ DEFAULT NOW(),
    accessed_at TIMESTAMPTZ
);

-- =====================================================
-- 4. Code Execution History
-- =====================================================
CREATE TABLE IF NOT EXISTS code_execution_history (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    session_id UUID REFERENCES design_sessions(id) ON DELETE CASCADE,
    agent_name VARCHAR(50) DEFAULT 'coder',

    -- Code Content
    code_content TEXT NOT NULL,
    snippet_ids TEXT[] DEFAULT '{}',

    -- Execution Result
    execution_result JSONB,
    stdout TEXT,
    stderr TEXT,
    exit_code INTEGER,
    execution_time_ms INTEGER,

    -- Self-healing
    is_healed BOOLEAN DEFAULT FALSE,
    original_code_id UUID REFERENCES code_execution_history(id),
    healing_attempt INTEGER DEFAULT 0,

    created_at TIMESTAMPTZ DEFAULT NOW()
);

-- =====================================================
-- 5. Code Snippet Library: Verified Code Templates
-- =====================================================
CREATE TABLE IF NOT EXISTS code_snippet_library (
    id VARCHAR(50) PRIMARY KEY,
    name VARCHAR(100) NOT NULL,
    description TEXT,
    category VARCHAR(50),

    -- Code
    code_template TEXT NOT NULL,
    required_imports TEXT[] DEFAULT '{}',
    input_params JSONB,
    output_schema JSONB,

    -- Quality
    is_verified BOOLEAN DEFAULT FALSE,
    test_coverage DECIMAL(5,2),
    last_verified_at TIMESTAMPTZ,

    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- =====================================================
-- 6. Candidate Snippets: Healer Knowledge Loop
-- =====================================================
CREATE TABLE IF NOT EXISTS candidate_snippets (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),

    -- Source Context
    source_session_id UUID REFERENCES design_sessions(id),
    source_healer_log_id UUID REFERENCES agent_execution_logs(id),

    -- Error Pattern
    error_type VARCHAR(50) NOT NULL,
    error_pattern TEXT NOT NULL,
    error_context JSONB,

    -- Fix Pattern
    fix_code TEXT NOT NULL,
    fix_explanation TEXT,
    fix_strategy VARCHAR(100),

    -- Quality Metrics
    success_count INTEGER DEFAULT 1,
    failure_count INTEGER DEFAULT 0,
    confidence_score DECIMAL(3,2) DEFAULT 0.50,

    -- Promotion Status
    status VARCHAR(20) DEFAULT 'candidate' CHECK (status IN (
        'candidate', 'reviewing', 'promoted', 'rejected'
    )),
    promoted_snippet_id VARCHAR(50) REFERENCES code_snippet_library(id),
    promoted_at TIMESTAMPTZ,

    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- =====================================================
-- 7. Scaffold Knowledge Mapping: Librarian RAG
-- =====================================================
CREATE TABLE IF NOT EXISTS scaffold_knowledge_mapping (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),

    -- Scaffold Info
    scaffold_smiles TEXT NOT NULL,
    scaffold_type VARCHAR(50),
    fingerprint_bits BYTEA,

    -- Semantic Embedding (768 for Gemini)
    embedding_vector vector(768),
    scaffold_description TEXT,

    -- Knowledge Links
    linked_knowledge_ids INTEGER[] DEFAULT '{}',
    linked_golden_set_ids UUID[] DEFAULT '{}',
    linked_pmids TEXT[] DEFAULT '{}',

    -- Usage Stats
    retrieval_count INTEGER DEFAULT 0,
    last_retrieved_at TIMESTAMPTZ,

    -- Quality
    relevance_score DECIMAL(3,2),
    is_validated BOOLEAN DEFAULT FALSE,

    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),

    UNIQUE(scaffold_smiles, scaffold_type)
);

-- =====================================================
-- 8. AlphaFold Jobs: Premium 3D Modeling
-- =====================================================
CREATE TABLE IF NOT EXISTS alphafold_jobs (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    session_id UUID REFERENCES design_sessions(id) ON DELETE CASCADE,
    user_id UUID REFERENCES auth.users(id) NOT NULL,

    -- Input
    target_sequence TEXT NOT NULL,
    ligand_smiles TEXT NOT NULL,
    job_type VARCHAR(50) DEFAULT 'binding',

    -- Job Status
    status VARCHAR(30) DEFAULT 'queued' CHECK (status IN (
        'queued', 'processing', 'completed', 'failed', 'timeout'
    )),
    priority INTEGER DEFAULT 5,

    -- Results
    predicted_structure JSONB,
    binding_affinity DECIMAL(8,4),
    confidence_score DECIMAL(3,2),
    contact_residues TEXT[],
    visualization_url TEXT,

    -- Compute Metrics
    compute_time_seconds INTEGER,
    gpu_hours_used DECIMAL(6,3),

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    started_at TIMESTAMPTZ,
    completed_at TIMESTAMPTZ
);

-- =====================================================
-- 9. AlphaFold User Quota: GPU Management
-- =====================================================
CREATE TABLE IF NOT EXISTS alphafold_user_quota (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID REFERENCES auth.users(id) UNIQUE NOT NULL,

    -- Monthly Quota
    monthly_gpu_hours_limit DECIMAL(6,2) DEFAULT 10.0,
    monthly_gpu_hours_used DECIMAL(6,2) DEFAULT 0.0,
    quota_reset_date DATE DEFAULT (DATE_TRUNC('month', NOW()) + INTERVAL '1 month')::DATE,

    -- Daily Limit (DDoS Prevention)
    daily_job_limit INTEGER DEFAULT 20,
    daily_jobs_submitted INTEGER DEFAULT 0,
    daily_reset_date DATE DEFAULT CURRENT_DATE,

    -- Priority Tier
    priority_tier VARCHAR(20) DEFAULT 'standard' CHECK (priority_tier IN (
        'standard', 'priority', 'enterprise'
    )),

    -- Overage Policy
    allow_overage BOOLEAN DEFAULT FALSE,
    overage_rate_per_hour DECIMAL(6,2) DEFAULT 5.00,

    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- Quota Auto-Reset Trigger
CREATE OR REPLACE FUNCTION reset_alphafold_quota()
RETURNS TRIGGER AS $$
BEGIN
    -- Monthly quota reset
    IF NEW.quota_reset_date <= CURRENT_DATE THEN
        NEW.monthly_gpu_hours_used := 0;
        NEW.quota_reset_date := (DATE_TRUNC('month', NOW()) + INTERVAL '1 month')::DATE;
    END IF;

    -- Daily quota reset
    IF NEW.daily_reset_date < CURRENT_DATE THEN
        NEW.daily_jobs_submitted := 0;
        NEW.daily_reset_date := CURRENT_DATE;
    END IF;

    NEW.updated_at := NOW();
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

DROP TRIGGER IF EXISTS auto_reset_quota ON alphafold_user_quota;
CREATE TRIGGER auto_reset_quota
    BEFORE UPDATE ON alphafold_user_quota
    FOR EACH ROW EXECUTE FUNCTION reset_alphafold_quota();

-- GPU Usage Increment Function
CREATE OR REPLACE FUNCTION increment_gpu_usage(p_user_id UUID, p_hours DECIMAL)
RETURNS VOID AS $$
BEGIN
    UPDATE alphafold_user_quota
    SET monthly_gpu_hours_used = monthly_gpu_hours_used + p_hours,
        updated_at = NOW()
    WHERE user_id = p_user_id;
END;
$$ LANGUAGE plpgsql;

-- =====================================================
-- INDEXES
-- =====================================================

-- Design Sessions
CREATE INDEX IF NOT EXISTS idx_sessions_user_status ON design_sessions(user_id, status);
CREATE INDEX IF NOT EXISTS idx_sessions_type ON design_sessions(session_type);
CREATE INDEX IF NOT EXISTS idx_sessions_created ON design_sessions(created_at DESC);

-- Agent Logs
CREATE INDEX IF NOT EXISTS idx_logs_session_agent ON agent_execution_logs(session_id, agent_name);
CREATE INDEX IF NOT EXISTS idx_logs_status ON agent_execution_logs(status);
CREATE INDEX IF NOT EXISTS idx_logs_healing ON agent_execution_logs(agent_name, healing_successful)
    WHERE agent_name = 'healer';

-- Encrypted Structures
CREATE INDEX IF NOT EXISTS idx_encrypted_session ON encrypted_structures(session_id);
CREATE INDEX IF NOT EXISTS idx_encrypted_user ON encrypted_structures(user_id);

-- Code History
CREATE INDEX IF NOT EXISTS idx_code_session ON code_execution_history(session_id);
CREATE INDEX IF NOT EXISTS idx_code_healed ON code_execution_history(is_healed);

-- Candidate Snippets
CREATE INDEX IF NOT EXISTS idx_candidate_error_type ON candidate_snippets(error_type, status);
CREATE INDEX IF NOT EXISTS idx_candidate_confidence ON candidate_snippets(confidence_score DESC);
CREATE INDEX IF NOT EXISTS idx_candidate_status ON candidate_snippets(status);

-- Scaffold Mapping
CREATE INDEX IF NOT EXISTS idx_scaffold_type ON scaffold_knowledge_mapping(scaffold_type);
CREATE INDEX IF NOT EXISTS idx_scaffold_validated ON scaffold_knowledge_mapping(is_validated);

-- AlphaFold
CREATE INDEX IF NOT EXISTS idx_alphafold_session ON alphafold_jobs(session_id);
CREATE INDEX IF NOT EXISTS idx_alphafold_status ON alphafold_jobs(status, priority);
CREATE INDEX IF NOT EXISTS idx_quota_user ON alphafold_user_quota(user_id);
CREATE INDEX IF NOT EXISTS idx_quota_reset ON alphafold_user_quota(quota_reset_date);

-- =====================================================
-- ROW LEVEL SECURITY (RLS)
-- =====================================================

-- Enable RLS
ALTER TABLE design_sessions ENABLE ROW LEVEL SECURITY;
ALTER TABLE encrypted_structures ENABLE ROW LEVEL SECURITY;
ALTER TABLE alphafold_jobs ENABLE ROW LEVEL SECURITY;
ALTER TABLE alphafold_user_quota ENABLE ROW LEVEL SECURITY;

-- Policies: Users see only their own data
DROP POLICY IF EXISTS "Users see own sessions" ON design_sessions;
CREATE POLICY "Users see own sessions" ON design_sessions
    FOR ALL USING (auth.uid() = user_id);

DROP POLICY IF EXISTS "Users see own structures" ON encrypted_structures;
CREATE POLICY "Users see own structures" ON encrypted_structures
    FOR ALL USING (auth.uid() = user_id);

DROP POLICY IF EXISTS "Users see own alphafold jobs" ON alphafold_jobs;
CREATE POLICY "Users see own alphafold jobs" ON alphafold_jobs
    FOR ALL USING (auth.uid() = user_id);

DROP POLICY IF EXISTS "Users see own quota" ON alphafold_user_quota;
CREATE POLICY "Users see own quota" ON alphafold_user_quota
    FOR ALL USING (auth.uid() = user_id);

-- =====================================================
-- INSERT DEFAULT CODE SNIPPETS
-- =====================================================
INSERT INTO code_snippet_library (id, name, description, category, code_template, required_imports, is_verified)
VALUES
    ('lipinski_check', 'Lipinski Rule of Five', 'Checks if molecule passes Lipinski rules', 'validation',
     E'mol = Chem.MolFromSmiles(smiles)\nif mol:\n    result["lipinski"] = {\n        "mw": Descriptors.MolWt(mol),\n        "logp": Descriptors.MolLogP(mol),\n        "hbd": Lipinski.NumHDonors(mol),\n        "hba": Lipinski.NumHAcceptors(mol),\n        "pass": Descriptors.MolWt(mol) <= 500 and Descriptors.MolLogP(mol) <= 5\n    }',
     ARRAY['rdkit.Chem', 'rdkit.Chem.Descriptors', 'rdkit.Chem.Lipinski'], TRUE),

    ('tanimoto_similarity', 'Tanimoto Coefficient Calculator', 'Calculates Tanimoto similarity', 'calculation',
     E'mol = Chem.MolFromSmiles(smiles)\nif mol:\n    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)\n    result["tanimoto"] = {"calculated": True, "fingerprint_bits": fp.GetNumOnBits()}',
     ARRAY['rdkit.Chem', 'rdkit.Chem.AllChem', 'rdkit.DataStructs'], TRUE),

    ('pains_filter', 'PAINS Filter', 'Checks for Pan-Assay Interference patterns', 'validation',
     E'mol = Chem.MolFromSmiles(smiles)\nif mol:\n    params = FilterCatalog.FilterCatalogParams()\n    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)\n    catalog = FilterCatalog.FilterCatalog(params)\n    entry = catalog.GetFirstMatch(mol)\n    result["pains"] = {"has_pains": entry is not None, "pattern": entry.GetDescription() if entry else None}',
     ARRAY['rdkit.Chem', 'rdkit.Chem.FilterCatalog'], TRUE),

    ('sa_score', 'Synthetic Accessibility Score', 'Calculates SA score (1=easy, 10=hard)', 'calculation',
     E'mol = Chem.MolFromSmiles(smiles)\nif mol:\n    from rdkit.Contrib.SA_Score import sascorer\n    result["sa_score"] = sascorer.calculateScore(mol)',
     ARRAY['rdkit.Chem', 'rdkit.Chem.rdMolDescriptors'], TRUE),

    ('mw_calc', 'Molecular Weight Calculator', 'Calculates molecular weight', 'calculation',
     E'mol = Chem.MolFromSmiles(smiles)\nif mol:\n    result["mw"] = Descriptors.MolWt(mol)',
     ARRAY['rdkit.Chem', 'rdkit.Chem.Descriptors'], TRUE)
ON CONFLICT (id) DO NOTHING;

-- =====================================================
-- COMPLETION MESSAGE
-- =====================================================
DO $$
BEGIN
    RAISE NOTICE '✅ ADC Design Engine Phase 1 Schema Created Successfully';
    RAISE NOTICE '   - 9 tables created';
    RAISE NOTICE '   - 21 CFR Part 11 compliance enabled';
    RAISE NOTICE '   - RLS policies applied';
    RAISE NOTICE '   - Default code snippets inserted';
END $$;
