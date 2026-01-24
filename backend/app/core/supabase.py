"""
Supabase Client
"""
from supabase import create_client, Client
from app.core.config import settings

class MockSupabaseClient:
    def table(self, name):
        return self
    def select(self, *args, **kwargs):
        return self
    def insert(self, *args, **kwargs):
        return self
    def update(self, *args, **kwargs):
        return self
    def upsert(self, *args, **kwargs):
        return self
    def eq(self, *args, **kwargs):
        return self
    def order(self, *args, **kwargs):
        return self
    def limit(self, *args, **kwargs):
        return self
    def execute(self):
        print("Warning: Executing on MockSupabaseClient (No DB connection)")
        return type('obj', (object,), {'data': [], 'count': 0})
    def rpc(self, *args, **kwargs):
        return self

def get_supabase_client() -> Client:
    """Supabase 클라이언트 반환"""
    url = settings.SUPABASE_URL
    if url and url.endswith("/"):
        url = url[:-1]
    key = settings.SUPABASE_SERVICE_KEY
    
    if not url or not key:
        print("Warning: SUPABASE_URL or SUPABASE_SERVICE_KEY is missing. Using Mock Client.")
        return MockSupabaseClient()
    
    try:
        return create_client(url, key)
    except Exception as e:
        print(f"Error creating Supabase client: {e}")
        return MockSupabaseClient()

supabase = get_supabase_client()

