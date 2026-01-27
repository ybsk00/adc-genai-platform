# Creative Biolabs Payload (Toxin) Crawler Implementation Plan

## 1. Overview
This crawler identifies and scrapes ADC Toxin payloads from Creative Biolabs and stores them in the `commercial_reagents` table. It is based on the architecture of `cb_antibody_crawler_lite.py`.

## 2. Target Source
- **Base URL:** `https://www.creative-biolabs.com/adc/classify-adc-toxins-6.htm`
- **Pagination:** Pages 1 to 8 (`?page=1` to `?page=8`)
- **Product Page Example:** `https://www.creative-biolabs.com/adc/auristatin-e-707.htm`

## 3. Data Extraction (Detail Page)
The crawler will extract data from the `<div class="proshow-unit">` section.

| Field | Source on Page | Database Column | Notes |
|---|---|---|---|
| **Product Number** | "Product No" or Title | `ambeed_cat_no` | **CRITICAL (Unique Key)**. Must be mapped (e.g., ADC-P-004). |
| **Product Name** | Page Title | `product_name` | Remove suffixes like " - Creative Biolabs" |
| **CAS NO** | `<li>CAS NO</li>` | `cas_number` | Important for future PubChem cleanup |
| **Classification** | `<li>Classification</li>` | `category` | Default to 'ADC Toxins' if missing. |
| **Molecular Weight** | `<li>Molecular Weight</li>` | `molecular_weight` | Store as text (e.g., "745.99") |
| **Purity** | `<li>Purity</li>` | `properties->>'purity'` | Store in JSONB properties |
| **Description** | `<li>Description</li>` | `properties->>'description'` | Store in JSONB properties |
| **Source** | Fixed Value | `source_name` | Always **'Creative Biolabs'** |

## 4. Database Schema (`commercial_reagents`)
The crawler will use the following mapping for the `commercial_reagents` table:

```sql
INSERT INTO commercial_reagents (
    ambeed_cat_no,
    product_name,
    category,
    cas_number,
    molecular_weight,
    source_name,
    product_url,
    properties,
    crawled_at
) VALUES (
    :cat_no,
    :name,
    :classification, -- or 'ADC Toxins'
    :cas,
    :mw,
    'Creative Biolabs',
    :url,
    :json_properties,
    NOW()
)
ON CONFLICT (ambeed_cat_no) DO UPDATE SET
    product_name = EXCLUDED.product_name,
    cas_number = EXCLUDED.cas_number,
    properties = EXCLUDED.properties,
    crawled_at = NOW();
```

## 5. Implementation Logic (`backend/cb_payload_crawler.py`)

### Class `CB_PayloadCrawler`
- **Init:** Setup Supabase, Playwright, Proxies (same as Antibody crawler).
- **Run Loop:**
    - Iterate `page` from 1 to 8.
    - URL: `https://www.creative-biolabs.com/adc/classify-adc-toxins-6.htm?page={page}`
    - **List Parsing:**
        - Select all product links (likely `a` tags containing `/adc/` and ending in `.htm` but NOT containing `classify-adc`).
    - **Detail Parsing (`process_product`):**
        - Navigate to product URL.
        - Extract fields using `page.evaluate`.
        - **Validation:** Ensure `ambeed_cat_no` is found. If not, log error and skip.
        - **Upsert:** Save to `commercial_reagents`.

## 6. Execution
- Run via: `python backend/cb_payload_crawler.py`

## 7. Future Work (Post-Processing)
- **AI Structure Fix:** Separate script/agent to iterate `commercial_reagents` where `cas_number` is present but `smiles_code` is NULL.
- Fetch SMILES from PubChem via API.
- Validate with LLM.
