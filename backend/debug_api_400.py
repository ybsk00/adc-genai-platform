import aiohttp
import asyncio
import json

API_BASE_URL = "https://clinicaltrials.gov/api/v2/studies"

async def test_api():
    search_term = 'Antibody Drug Conjugate'
    last_update_date = "2026-01-25"
    status_filter = [
        "COMPLETED", "TERMINATED", "WITHDRAWN", 
        "SUSPENDED", "RECRUITING", "ACTIVE_NOT_RECRUITING",
        "ENROLLING_BY_INVITATION", "NOT_YET_RECRUITING", "UNKNOWN"
    ]
    
    # API v2 syntax for date filtering
    query_term = f"{search_term} AND AREA[LastUpdatePostDate]RANGE[{last_update_date},MAX]"
    
    params = {
        "query.term": query_term,
        "filter.overallStatus": ",".join(status_filter),
        "pageSize": 10,
        "format": "json"
    }
    
    print(f"Testing API with params: {params}")
    
    async with aiohttp.ClientSession() as session:
        async with session.get(API_BASE_URL, params=params) as response:
            print(f"Status: {response.status}")
            print(f"Headers: {dict(response.headers)}")
            try:
                data = await response.json()
                print(f"Response Body: {json.dumps(data, indent=2)}")
            except:
                text = await response.text()
                print(f"Response Text: {text}")

if __name__ == "__main__":
    asyncio.run(test_api())
