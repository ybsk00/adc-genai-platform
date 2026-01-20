from app.core.supabase import supabase
import json

def check_db():
    tables = ["golden_set", "golden_set_library", "knowledge_base", "data_sync_logs"]
    results = {}
    
    for table in tables:
        try:
            res = supabase.table(table).select("count", count="exact").limit(1).execute()
            results[table] = {"exists": True, "count": res.count}
        except Exception as e:
            results[table] = {"exists": False, "error": str(e)}
            
    print(json.dumps(results, indent=2))

if __name__ == "__main__":
    check_db()
