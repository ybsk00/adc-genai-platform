from app.core.supabase import supabase
import sys

def add_column():
    try:
        # Supabase client doesn't have a direct SQL execution method unless RPC is set up.
        # Let's try to see if we can use a common RPC name like 'exec_sql' if it exists.
        # If not, we might need to rely on the user to run it or find another way.
        
        sql = "ALTER TABLE sync_jobs ADD COLUMN IF NOT EXISTS last_processed_page INTEGER DEFAULT 0;"
        
        # Try a common RPC for executing SQL (if the user has set it up)
        try:
            res = supabase.rpc('exec_sql', {'sql_query': sql}).execute()
            print("✅ Successfully added column via RPC 'exec_sql'")
            return
        except Exception as e:
            print(f"⚠️ RPC 'exec_sql' failed: {e}")

        try:
            res = supabase.rpc('run_sql', {'sql': sql}).execute()
            print("✅ Successfully added column via RPC 'run_sql'")
            return
        except Exception as e:
            print(f"⚠️ RPC 'run_sql' failed: {e}")

        print("❌ Could not add column via standard RPCs. Please add 'last_processed_page' (INTEGER) to 'sync_jobs' table manually.")
        
    except Exception as e:
        print(f"🔥 Error: {e}")

if __name__ == "__main__":
    add_column()
