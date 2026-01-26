# ADC 플랫폼 데이터 수집 가속화 및 고도화 계획 (v2)

## 1. 구현 완료 내역 (Implementation Status)

### A. 연도별 분할 수집 (Yearly Chunking Strategy)
- **기존:** 최근 30일 데이터만 단일 호출.
- **변경:** 2026년부터 2010년까지 역순으로 **연도별 루프(Year Loop)**를 돌며 수집.
- **효과:** Entrez API의 1회 호출 제한(10,000건)을 우회하여 이론상 무제한 수집 가능.

### B. Mega-Net 쿼리 확장 (Query Expansion)
- **기존:** 단순 MeSH Term 및 ADC 키워드.
- **변경:** 
    - **Suffix Wildcards:** `*-mab`, `*-tin`, `*-can`, `*-tecan`, `*-vedotin`, `*-deruxtecan` 등 포괄적 접미사 추가.
    - **MeSH Terms:** `Immunoconjugates`, `Bispecific Antibodies`, `Checkpoint Inhibitor` 등 관련성 높은 상위 개념 포함.
    - **검색 범위:** 기존 대비 약 5배 이상의 검색 커버리지 확보.

### C. 데이터 지능화 (Inteligent Extraction)
- **4단 분리:** AI 프롬프트를 강화하여 ADC 구조를 **Target / Antibody / Linker / Payload**로 명확히 분리.
- **임상 수치 정형화:** ORR(%), PFS(mo), OS(mo) 수치를 텍스트가 아닌 **숫자(Float)** 필드로 추출.
- **NCT ID 강제 주입:** AI가 놓칠 수 있는 임상시험 번호를 **Regex(정규식)**로 우선 추출하여 강제로 `properties`에 주입.

### D. 중복 관리 (Upsert Logic)
- **기존:** PMID 중복 시 단순히 Skip (정보 갱신 불가).
- **변경:** 중복 발견 시 기존 레코드의 `properties` 및 `analysis` 필드를 최신 정보로 **Update (Upsert)**.

---

## 2. 추가 모수 확보를 위한 제언 (Future Expansion Ideas)

PubMed 외에 데이터 모수를 획기적으로 늘리기 위한 추가 전략입니다.

### 1. ClinicalTrials.gov API 직접 연동
- **설명:** 논문에서 NCT ID를 추출하는 것을 넘어, ClinicalTrials.gov API를 직접 찔러서 "ADC" 관련 모든 임상시험을 가져옵니다.
- **장점:** 논문이 아직 나오지 않은 **진행 중(Recruiting)**인 임상 정보를 가장 빠르게 확보 가능.
- **예상 모수:** 약 2,000~3,000건 추가 가능.

### 2. 주요 학회 초록 (Conference Abstracts) 크롤링
- **타겟:** ASCO (미국임상종양학회), ESMO (유럽종양학회), AACR (미국암학회).
- **이유:** 최신 ADC 임상 결과는 정식 논문(PubMed) 등재 전, 학회에서 먼저 초록 형태로 발표됩니다 (6개월~1년 선행).
- **방법:** 각 학회 아카이브 사이트에 대한 Custom Python Crawler 개발 필요.

### 3. 특허 데이터베이스 (Patent Mining)
- **타겟:** Google Patents Public Data (BigQuery) 또는 USPTO/WIPO API.
- **이유:** 초기 단계의 새로운 Payload 기술이나 Linker 플랫폼은 특허에 먼저 등장합니다.
- **검색어:** "Antibody-Drug Conjugate" AND ("Linker" OR "Payload").

### 4. 제약사 파이프라인 페이지 모니터링
- **타겟:** Daiichi Sankyo, Seagen (Pfizer), Roche, AbbVie, Gilead 등 Top-tier ADC 플레이어의 공식 'Pipeline' 웹페이지.
- **방법:** `Playwright` 등을 이용한 주기적 스크래핑. 개발 코드명(예: DS-xxxx)과 타겟 정보를 매핑하는 데 매우 유용.
