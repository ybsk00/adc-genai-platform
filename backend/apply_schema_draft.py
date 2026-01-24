import os
import logging
from dotenv import load_dotenv
from app.core.supabase import supabase

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Load env
load_dotenv(os.path.join(os.path.dirname(__file__), '.env'))

def apply_schema():
    logger.info("Applying schema updates...")
    
    # Read SQL file
    with open('backend/update_commercial_reagents_schema.sql', 'r', encoding='utf-8') as f:
        sql = f.read()
        
    # Split by statement if needed, or try running as one block if supported by postgres rpc or direct sql function
    # Supabase-py doesn't have a direct 'query' method for raw SQL unless we use a stored procedure or pg driver.
    # However, we can use the `rpc` method if we have a function to exec sql, OR we can try to use a direct connection if we had psycopg2.
    # But wait, the user instructions say "SQL을 통한 차원 검증" using `SELECT ...`.
    # If I can't run psql, I might need to use a python script with psycopg2 if installed.
    # Let's check if psycopg2 is in requirements or used in the project.
    # I see `app.core.supabase` which uses `supabase` client (postgrest). Postgrest doesn't support raw SQL execution easily unless enabled.
    
    # Alternative: Use the `postgres` connection string from .env and `psycopg2` or `asyncpg`.
    # Let's try to find a way to execute raw SQL.
    # If `psycopg2` is not available, I might have to rely on the user to run it or use a pre-existing `exec_sql` function in the DB if it exists.
    
    # Let's check if there is a `exec_sql` RPC function.
    try:
        # Try to call a hypothetical exec_sql function
        # supabase.rpc('exec_sql', {'query': sql}).execute()
        # If that doesn't exist, we might be stuck without psql.
        pass
    except Exception:
        pass

    # Actually, the best way if psql failed is to use python with psycopg2 if available.
    try:
        import psycopg2
        db_url = os.getenv("DATABASE_URL") # Check if this env var exists
        if not db_url:
            # Try to construct from other vars or finding it in .env
            # The user has a .env file.
            pass
            
        conn = psycopg2.connect(db_url)
        cur = conn.cursor()
        cur.execute(sql)
        conn.commit()
        cur.close()
        conn.close()
        logger.info("✅ Schema applied successfully via psycopg2.")
        return
    except ImportError:
        logger.error("❌ psycopg2 not installed.")
    except Exception as e:
        logger.error(f"❌ Failed to apply schema: {e}")

if __name__ == "__main__":
    apply_schema()
