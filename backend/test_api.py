"""
API 테스트 스크립트
"""
import aiohttp
import asyncio

async def test_clinical_trials():
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
        'Accept': 'application/json'
    }
    async with aiohttp.ClientSession(headers=headers) as session:
        async with session.get(
            'https://clinicaltrials.gov/api/v2/studies',
            params={'query.term': 'ADC', 'pageSize': 5, 'format': 'json'}
        ) as res:
            print('Status:', res.status)
            if res.status == 200:
                data = await res.json()
                print('Studies found:', len(data.get('studies', [])))
            else:
                text = await res.text()
                print('Error:', text[:500])

asyncio.run(test_clinical_trials())
