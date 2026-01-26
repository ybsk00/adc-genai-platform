from app.core.supabase import supabase
import asyncio

async def check_columns():
    try:
        # Try to select the new columns
        res = supabase.table("golden_set_library").select("dar, patient_count, adverse_events_grade3_pct, target_symbol").limit(1).execute()
        print("✅ Columns exist.")
    except Exception as e:
        print(f"❌ Error selecting columns: {e}")

if __name__ == "__main__":
    asyncio.run(check_columns())
