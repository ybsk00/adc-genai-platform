from app.core.supabase import supabase
import json

try:
    # Get one record to see the columns
    res = supabase.table('commercial_reagents').select('*').limit(1).execute()
    if res.data:
        print("Columns found in commercial_reagents:")
        print(", ".join(list(res.data[0].keys())))
    else:
        print("No records found in commercial_reagents to infer schema.")
except Exception as e:
    print(f"Error: {e}")
