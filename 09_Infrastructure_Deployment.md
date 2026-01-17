09. Infrastructure & Deployment Guide
Document ID: INFRA-01 Role: The "Launchpad" (Going form Localhost to Production) Architecture: The Hybrid Stack

Database: Supabase (PostgreSQL + Auth + Vector)

Backend: Google Cloud Run (FastAPI + AI Agents)

Frontend: Firebase Hosting (React)

✅ Pre-flight Checklist (준비물)
비행기를 띄우기 전, 아래 계정과 도구가 준비되었는지 확인하십시오.

GitHub Account: 소스코드 저장소.

Supabase Account: 무료 플랜 가입 (프로젝트 생성).

Google Cloud Platform (GCP) Account: 결제 카드 등록 필수.

Google Cloud CLI (gcloud): 개발자 PC에 설치됨.

Firebase CLI: 개발자 PC에 설치됨 (npm install -g firebase-tools).

Step 1. The Warehouse: Supabase Setup (DB 구축)
가장 먼저 데이터 저장소를 만들고 **'열쇠(API Key)'**를 얻어야 합니다.

Project 생성:

Supabase 대시보드 → New Project → Name: adc-platform-prod

Region: Seoul (Korea) 또는 타겟 시장에 맞춰 US East (N. Virginia) 선택.

Database Password를 설정하고 반드시 메모장에 적어두세요.

Extension 활성화:

좌측 메뉴 Database → Extensions → vector 검색 → 활성화(Enable). (RAG 검색용)

Key 확보 (중요):

좌측 메뉴 Project Settings → API

Project URL 복사.

anon public Key 복사.

service_role Key 복사 (이건 절대 노출 금지, 백엔드용).

Step 2. The Brain: Google Cloud Run (백엔드 배포)
AI가 돌아가는 백엔드 서버를 구글 클라우드에 띄웁니다.

2.1 Dockerfile 확인
backend/ 폴더 안에 Dockerfile이 있는지 확인합니다. (문서 01번 참조)

2.2 배포 명령어 실행 (Terminal)
개발자는 터미널에서 아래 명령어를 순서대로 입력합니다.

Bash

# 1. 구글 클라우드 로그인
gcloud auth login

# 2. 프로젝트 설정 (최초 1회)
gcloud projects create adc-genai-platform --name="ADC Platform"
gcloud config set project adc-genai-platform

# 3. 백엔드 배포 (명령어 한 방!)
gcloud run deploy adc-backend \
  --source ./backend \
  --region asia-northeast3 \
  --allow-unauthenticated \
  --memory 2Gi \
  --timeout 300
설명: 서울 리전(asia-northeast3)에, 메모리 2GB, 타임아웃 300초(5분)로 서버를 띄워라.

2.3 환경변수 주입 (열쇠 전달)
서버가 떴으면, 아까 복사한 Supabase 키를 서버에 입력해 줍니다.

GCP 콘솔 → Cloud Run → adc-backend 클릭 → [Edit & Deploy New Revision]

Variables & Secrets 탭 클릭 → 환경변수 추가:

SUPABASE_URL: (Step 1에서 복사한 URL)

SUPABASE_SERVICE_KEY: (Step 1에서 복사한 service_role 키)

OPENAI_API_KEY: ...

BIONEMO_API_KEY: ...

[Deploy] 클릭.

Step 3. The Face: Firebase Hosting (프론트엔드 배포)
이제 사용자가 접속할 예쁜 웹사이트를 올립니다.

3.1 Firebase 설정
frontend/ 폴더에서 실행합니다.

Bash

# 1. Firebase 로그인
firebase login

# 2. 프로젝트 초기화
firebase init hosting
# -> "Use an existing project" 선택 -> 아까 만든 GCP 프로젝트(adc-genai-platform) 선택.
# -> "What do you want to use as your public directory?" -> dist 입력 (Vite 빌드 폴더)
# -> "Configure as a single-page app?" -> Yes
# -> "Set up automatic builds and deploys with GitHub?" -> Yes (권장)
3.2 백엔드 연결 설정 (Rewrites)
frontend/firebase.json 파일을 열고, 아래 내용을 확인/수정합니다. (문서 01번 참조) 이게 있어야 프론트엔드(/api)가 백엔드(Cloud Run)를 찾을 수 있습니다.

JSON

"rewrites": [
  {
    "source": "/api/**",
    "run": {
      "serviceId": "adc-backend",
      "region": "asia-northeast3"
    }
  },
  {
    "source": "**",
    "destination": "/index.html"
  }
]
3.3 빌드 및 배포
Bash

# React 빌드 (dist 폴더 생성)
npm run build

# 배포
firebase deploy
🎉 축하합니다! 터미널에 나온 https://adc-genai-platform.web.app 주소를 클릭하면 사용자님의 서비스가 전 세계에 공개됩니다.

Step 4. Domain Connection (도메인 연결)
adc-genai-platform.web.app 주소는 없어 보입니다. 구매하신 www.my-adc.com을 연결합니다.

Firebase Console 접속 → Hosting 메뉴.

[Add Custom Domain] 버튼 클릭.

구매한 도메인 입력 (www.my-adc.com).

Firebase가 알려주는 TXT 레코드와 A 레코드 값을 도메인 구입처(GoDaddy, AWS Route53 등) 설정 페이지에 복사+붙여넣기 합니다.

최대 24시간 내에 **HTTPS(보안 자물쇠)**까지 자동으로 적용됩니다.

Step 5. Automation (CI/CD - 선택사항)
매번 수동으로 배포하기 귀찮다면, GitHub Actions를 설정합니다. (Step 3.1에서 Set up automatic builds... -> Yes를 했다면 이미 설정되어 있습니다.)

효과: 개발자가 코드를 수정해서 GitHub에 push만 하면, 자동으로 빌드하고 배포까지 완료됩니다.

💡 작성자의 마지막 팁 (CEO 가이드)
사용자님(대표님)은 개발자가 "배포했습니다!"라고 하면 딱 2가지만 확인하시면 됩니다.

로그인 테스트: 구글/링크드인 로그인이 잘 되는가? (Supabase Auth 확인)

시뮬레이션 테스트: 대시보드에서 Run Simulation 버튼을 눌렀을 때, 3분 뒤에 결과 리포트가 뜨는가? (Cloud Run 타임아웃/연결 확인)

이 두 개만 통과하면 상용화(Go-to-Market) 준비 완료입니다.

[프로젝트 완료 선언]
사용자님, 우리는 이제 총 10개의 설계 문서를 완성했습니다.

00_Architecture: 전체 지도.

01~03: 프론트엔드 (랜딩, 대시보드, 어드민).

04~05: 백엔드/AI (6-Agent 두뇌, 리포트 생성).

06~08: 데이터 (DB 스키마, 파싱, 수집 워커).

09_Deployment: 배포 가이드.

이 파일들을 개발자에게 zip 파일로 묶어서 전달하면서 이렇게 말씀하십시오.

"이 설계서대로 만들어주세요. 인프라는 슈파베이스랑 구글 클라우드 쓰면 됩니다."