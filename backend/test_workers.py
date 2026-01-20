import asyncio
import httpx
import time

BASE_URL = "http://localhost:8000/api/scheduler"

async def test_sync():
    async with httpx.AsyncClient() as client:
        # 1. PubMed 동기화 트리거
        print("Triggering PubMed sync...")
        res = await client.post(f"{BASE_URL}/sync/pubmed")
        pubmed_job = res.json()["job_id"]
        print(f"PubMed Job ID: {pubmed_job}")

        # 2. OpenFDA 동기화 트리거
        print("Triggering OpenFDA sync...")
        res = await client.post(f"{BASE_URL}/sync/openfda")
        openfda_job = res.json()["job_id"]
        print(f"OpenFDA Job ID: {openfda_job}")

        # 3. 상태 모니터링 (404 오류 체크 포함)
        for _ in range(10):
            time.sleep(5)
            print(f"\nChecking status at {time.strftime('%H:%M:%S')}...")
            
            for job_id in [pubmed_job, openfda_job]:
                try:
                    status_res = await client.get(f"{BASE_URL}/sync/{job_id}")
                    status_data = status_res.json()
                    print(f"Job {job_id}: status={status_data['status']}, found={status_data['records_found']}, drafted={status_data['records_drafted']}")
                except Exception as e:
                    print(f"Error checking job {job_id}: {e}")

if __name__ == "__main__":
    try:
        asyncio.run(test_sync())
    except Exception as e:
        print(f"Test failed: {e}. (Make sure the backend server is running at localhost:8000)")
