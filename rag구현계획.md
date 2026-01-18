# [작업지시서] ADC 플랫폼 하이브리드 데이터 파이프라인 & RAG 구축 가이드

**문서 번호:** SPEC-DATA-02  
**작성 일자:** 2026-01-18  
**수신:** 개발팀  
**목표:** 고품질 ADC(Antibody-Drug Conjugate) 연구 데이터를 확보하기 위한 **수동 업로드** 및 **자동 수집(크롤링)** 이원화 시스템 구축

---

## 1. 시스템 개요 (System Architecture)

본 시스템은 **RAG(검색 증강 생성)**의 정확도를 극대화하기 위해 두 가지 경로로 데이터를 수집하고, **Vector DB(Supabase)**로 통합 관리한다.

1.  **Manual Pipeline (Admin Upload):** 관리자(CEO)가 보유한 유료 논문, 내부 실험 데이터(Excel/PDF)를 직접 업로드.
2.  **Auto Pipeline (Intelligent Crawler):** Python 봇이 매일 Open Access 논문을 수집하고, **AI(Gemini 2.0)**가 1차 검수 후 적재.

---

## 2. 상세 구현 요건 (Detailed Requirements)

### 🟢 A. 수동 업로드 파이프라인 (Manual Upload)
*관리자가 고품질의 핵심 데이터(Core Data)를 직접 주입하는 경로.*

#### 1. Admin UI 개발
* **파일 업로드 페이지:** Drag & Drop 지원.
* **지원 형식:** PDF (논문), Excel/CSV (실험 데이터).
* **메타데이터 입력:** 파일 업로드 시 `Title`, `Author`, `Category(특허/논문/실험)`, `Security Level(공개/비밀)` 태그 입력 기능.

#### 2. 처리 로직 (Processor)
* **Trigger:** 파일 업로드 즉시 실행.
* **Text Extraction:** PDF 내 텍스트 및 **표(Table)** 데이터 추출.
* **Embedding:** 추출된 텍스트를 청크(Chunk) 단위로 쪼개어 임베딩 후 Supabase에 저장.
* **Model:** 속도와 긴 문맥 처리를 위해 **Gemini 2.0 Flash** 사용.

---

### 🔵 B. 자동 수집 파이프라인 (Auto Crawling)
*양적 확대를 위해 외부 무료 논문을 자동으로 수집하는 경로. (노이즈 필터링 필수)*

#### 1. Crawler (Python Script)
* **Target:** PubMed Central (PMC), bioRxiv (Open Access Only).
* **Schedule:** 매일 자정 (00:00 KST) 실행 via Cloud Scheduler.
* **1차 필터 (Query Tuning):**
    * 단순 `ADC` 검색 금지 (노이즈 발생).
    * **필수 쿼리:** `("Antibody-Drug Conjugate" OR "ADC") AND ("Linker" OR "Payload" OR "Conjugation")`
    * **제외 쿼리:** `NOT ("Monotherapy" AND "Trastuzumab")` (단독 항체 치료 제외)

#### 2. AI Gatekeeper (2차 지능형 필터) ⭐ **(핵심)**
* **문제:** 검색어만으로는 '일반 항암제'나 '단순 항체 치료' 논문이 섞여 들어옴.
* **해결:** 수집된 논문의 초록(Abstract)을 **Gemini 2.0 Flash**에게 먼저 보내 검사함.
* **Prompt Logic:**
    > "이 논문이 항체(Antibody)와 약물(Payload)이 링커(Linker)로 결합된 기술을 다루고 있는가? 단순 항체 치료제라면 False를 반환하라."
* **Action:** `True` 판정을 받은 논문만 PDF 다운로드 및 DB 적재. `False`는 폐기(Log 남김).

---

## 3. 통합 RAG 구조 (Integrated RAG)

데이터 소스가 어디든(Manual/Auto), RAG는 동일한 Vector DB를 참조하여 답변을 생성한다.

1.  **Retriever (검색기):**
    * 사용자 질문 시 Supabase에서 관련 청크(Chunk) 검색.
    * *예: "링커 기술 트렌드 알려줘" -> Manual로 올린 특허 + Auto로 긁어온 논문 모두 검색.*

2.  **Generator (답변 생성기):**
    * **Role:** 검색된 자료들을 종합하여 최종 리포트 작성.
    * **Model:** **GPT-4o** (복잡한 추론 및 보고서 작성용).
    * **Process:** [사용자 질문] + [검색된 논문 요약] + [검색된 내부 실험 데이터] -> 최종 답변 생성.

---

## 4. 기술 스택 및 환경 설정 (.env)

* **Language:** Python 3.9+
* **Database:** Supabase (pgvector)
* **AI Models:**
    * **Preprocessing / Gatekeeper:** Google Gemini 2.0 Flash (`gemini-2.0-flash-exp`)
    * **Final Reporting:** OpenAI GPT-4o (`gpt-4o`)
* **Crawler Lib:** `Biopython`, `BeautifulSoup4`, `Selenium`

---

## 5. 개발자 체크리스트

1. [ ] **Admin 페이지:** 내가 가진 PDF 파일을 업로드하면 DB에 잘 들어가는가?
2. [ ] **크롤러 필터:** '트라스주맙' 단독 치료 논문이 자동으로 걸러지는가? (AI Gatekeeper 작동 확인)
3. [ ] **통합 검색:** RAG 질의 시, 내가 올린 파일과 크롤링한 파일이 동시에 검색 결과로 나오는가?