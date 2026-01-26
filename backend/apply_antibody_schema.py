import os
import logging
from dotenv import load_dotenv
from supabase import create_client

load_dotenv()
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("ApplySchema")

def apply_schema():
    url = os.getenv("SUPABASE_URL")
    key = os.getenv("SUPABASE_SERVICE_KEY")
    
    if not url or not key:
        logger.error("Missing Supabase credentials.")
        return

    # Read SQL file
    try:
        with open("create_antibody_tables.sql", "r", encoding="utf-8") as f:
            sql_content = f.read()
    except FileNotFoundError:
        logger.error("SQL file not found.")
        return

    # Split commands (simple split by ;)
    # Note: Supabase-py doesn't support raw SQL execution directly via client in all versions.
    # We will try to use the `rpc` function if a generic sql exec function exists, 
    # OR we use `postgrest` client directly if possible. 
    # HOWEVER, standard supabase-py is for data manipulation. 
    # FOR DDL (CREATE TABLE), we usually need the SQL Editor or a postgres connection (psycopg2).
    # Assuming 'postgres' connection is NOT available here (only REST API).
    
    # WORKAROUND: We will try to use a pre-existing RPC function 'exec_sql' if it exists (common pattern).
    # If not, we will inform the user to run the SQL in the Dashboard.
    
    client = create_client(url, key)
    
    try:
        # Attempt to run via RPC (if 'exec_sql' or similar exists)
        # Often projects have a helper for this. If not, we might fail.
        logger.info("Attempting to apply schema via RPC 'exec_sql'...")
        res = client.rpc("exec_sql", {"sql": sql_content}).execute()
        logger.info(f"Schema applied via RPC: {res}")
    except Exception as e:
        logger.warning(f"RPC 'exec_sql' failed or not found: {e}")
        logger.info("⚠️ Supabase REST API does not support CREATE TABLE directly.")
        logger.info("👉 Please copy the content of 'backend/create_antibody_tables.sql' and run it in the Supabase SQL Editor.")

if __name__ == "__main__":
    apply_schema()
