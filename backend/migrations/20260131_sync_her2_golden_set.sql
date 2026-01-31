-- ============================================================================
-- Migration: Sync HER2/Trastuzumab data to golden_set_library
-- Date: 2026-01-31
-- Purpose: Ensure golden_set_library has complete HER2 approved data
--          with linker_type, dar, orr_pct, pfs_months populated
-- ============================================================================

-- 1. Check current HER2 data in golden_set_library
-- SELECT name, status, target_1, linker_type, dar, orr_pct, pfs_months, os_months
-- FROM golden_set_library WHERE target_1 = 'HER2';

-- 2. Update existing HER2 records: fill missing linker/DAR fields
-- Trastuzumab emtansine (T-DM1 / Kadcyla): FDA approved 2013
UPDATE golden_set_library
SET
    linker_type = COALESCE(linker_type, 'non-cleavable'),
    dar = COALESCE(dar, 3.5),
    status = 'approved'
WHERE target_1 = 'HER2'
  AND name ILIKE '%trastuzumab emtansine%' OR name ILIKE '%kadcyla%' OR name ILIKE '%T-DM1%';

-- Trastuzumab deruxtecan (T-DXd / Enhertu): FDA approved 2019
UPDATE golden_set_library
SET
    linker_type = COALESCE(linker_type, 'cleavable'),
    dar = COALESCE(dar, 8),
    status = 'approved'
WHERE target_1 = 'HER2'
  AND (name ILIKE '%trastuzumab deruxtecan%' OR name ILIKE '%enhertu%' OR name ILIKE '%T-DXd%');

-- 3. Insert canonical HER2 ADCs if not already present
-- T-DM1 (Kadcyla)
INSERT INTO golden_set_library (
    name, target_1, category, status, linker_type, dar,
    orr_pct, pfs_months, os_months, outcome_type, description
)
SELECT
    'Trastuzumab emtansine (T-DM1)',
    'HER2',
    'ADC',
    'approved',
    'non-cleavable',
    3.5,
    43.6,
    9.6,
    30.9,
    'Approved',
    'FDA-approved ADC for HER2+ metastatic breast cancer. Maytansine derivative (DM1) payload with non-cleavable thioether linker.'
WHERE NOT EXISTS (
    SELECT 1 FROM golden_set_library
    WHERE target_1 = 'HER2'
      AND (name ILIKE '%T-DM1%' OR name ILIKE '%kadcyla%' OR name ILIKE '%trastuzumab emtansine%')
);

-- T-DXd (Enhertu)
INSERT INTO golden_set_library (
    name, target_1, category, status, linker_type, dar,
    orr_pct, pfs_months, os_months, outcome_type, description
)
SELECT
    'Trastuzumab deruxtecan (T-DXd)',
    'HER2',
    'ADC',
    'approved',
    'cleavable',
    8,
    79.7,
    25.1,
    NULL,
    'Approved',
    'FDA-approved ADC for HER2+ metastatic breast cancer. Topoisomerase I inhibitor (DXd) payload with cleavable tetrapeptide linker. DESTINY-Breast03 trial.'
WHERE NOT EXISTS (
    SELECT 1 FROM golden_set_library
    WHERE target_1 = 'HER2'
      AND (name ILIKE '%T-DXd%' OR name ILIKE '%enhertu%' OR name ILIKE '%trastuzumab deruxtecan%')
);

-- 4. Update Zongertinib (if exists) with linker/DAR info
-- Zongertinib is a TKI (small molecule), not ADC - but if it's in the DB:
UPDATE golden_set_library
SET
    linker_type = COALESCE(linker_type, 'N/A (TKI)'),
    dar = COALESCE(dar, NULL)
WHERE target_1 = 'HER2'
  AND name ILIKE '%zongertinib%';

-- 5. Verify results
-- SELECT name, status, target_1, linker_type, dar, orr_pct, pfs_months, os_months
-- FROM golden_set_library WHERE target_1 = 'HER2' ORDER BY orr_pct DESC NULLS LAST;
