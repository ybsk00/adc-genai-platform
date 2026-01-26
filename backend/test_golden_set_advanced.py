"""
Test Golden Set Advanced Extraction
- 이중항체 (Bispecifics)
- 독성 수치 (Grade 3+ AEs)
- DAR (Drug-to-Antibody Ratio)
- 신뢰도 점수
"""
import asyncio
import logging
from app.services.ai_refiner import ai_refiner

# Setup Logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def test_extraction():
    # 1. Bispecific Antibody Test Case
    bispecific_record = {
        "id": "test_bsab_001",
        "name": "Amivantamab",
        "enrichment_source": "clinical_trials",
        "properties": {
            "brief_summary": "Amivantamab is a bispecific antibody targeting EGFR and MET. It is being evaluated in patients with NSCLC.",
            "overall_status": "Active",
            "phase": "Phase 2"
        }
    }
    
    logger.info("🧪 Testing Bispecific Extraction...")
    res_bsab = await ai_refiner.refine_single_record(bispecific_record)
    logger.info(f"Result: {res_bsab}")
    
    assert res_bsab.get("target_1") in ["EGFR", "MET"], "Target 1 extraction failed"
    assert res_bsab.get("target_2") in ["EGFR", "MET"], "Target 2 extraction failed"
    assert res_bsab.get("antibody_format") in ["Bispecific", "bispecific"], "Format extraction failed"

    # 2. ADC Toxicity & DAR Test Case
    adc_record = {
        "id": "test_adc_001",
        "name": "Trastuzumab deruxtecan",
        "enrichment_source": "clinical_trials",
        "properties": {
            "brief_summary": "T-DXd is an ADC with a drug-to-antibody ratio (DAR) of approximately 8. In the DESTINY-Breast03 trial, Grade 3 or higher adverse events occurred in 45% of patients.",
            "overall_status": "Completed",
            "phase": "Phase 3"
        }
    }
    
    logger.info("\n🧪 Testing ADC Metrics Extraction...")
    res_adc = await ai_refiner.refine_single_record(adc_record)
    logger.info(f"Result: {res_adc}")
    
    dar_val = res_adc.get("dar")
    assert dar_val is not None and float(dar_val) == 8.0, f"DAR extraction failed: {dar_val}"
    # Note: LLM might extract 45 or 45% or 0.45. Adjust assertion based on prompt instruction (number only).
    # Prompt says "number only".
    assert res_adc.get("adverse_events_grade3_pct") is not None, "AE extraction failed"
    
    logger.info("\n✅ All Advanced Extraction Tests Passed!")

if __name__ == "__main__":
    asyncio.run(test_extraction())
