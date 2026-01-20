from app.core.supabase import supabase

def create_settings_table():
    try:
        # RPC를 사용하여 테이블 생성 시도 (Supabase SQL Editor 권장이나 API로 시도)
        # 직접적인 CREATE TABLE은 REST API로 불가능하므로, 
        # 사용자가 이미 schema.sql을 실행했다고 가정하거나, 
        # 여기서는 설정을 저장할 때 테이블이 없으면 에러가 나도록 하고 
        # 사용자에게 SQL 실행을 요청하는 것이 안전함.
        # 하지만 자동화를 위해 rpc('exec_sql') 등이 있다면 사용 가능.
        
        # 일단 테이블 존재 여부 확인
        res = supabase.table("data_source_settings").select("count", count="exact").limit(1).execute()
        print("Table 'data_source_settings' already exists.")
    except Exception as e:
        print(f"Table might not exist or error: {e}")
        print("Please run the following SQL in Supabase SQL Editor:")
        print("""
        CREATE TABLE IF NOT EXISTS data_source_settings (
            source_id TEXT PRIMARY KEY,
            auto_sync BOOLEAN DEFAULT FALSE,
            sync_interval_hours INTEGER DEFAULT 24,
            last_run_at TIMESTAMPTZ,
            next_run_at TIMESTAMPTZ,
            updated_at TIMESTAMPTZ DEFAULT NOW()
        );
        """)

if __name__ == "__main__":
    create_settings_table()
