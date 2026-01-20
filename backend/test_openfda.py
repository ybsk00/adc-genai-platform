"""
OpenFDA API 테스트 스크립트
"""
import httpx

# OpenFDA API 테스트
res = httpx.get(
    'https://api.fda.gov/drug/label.json',
    params={'search': 'openfda.generic_name:vedotin', 'limit': 5},
    timeout=30
)
print('Status:', res.status_code)
if res.status_code == 200:
    data = res.json()
    print('Results:', len(data.get('results', [])))
    for r in data.get('results', []):
        print(' -', r.get('openfda', {}).get('brand_name', ['N/A'])[0])
else:
    print('Error:', res.text[:500])
