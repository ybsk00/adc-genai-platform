01. Frontend Design: Landing, Auth & Payments
Document ID: FE-01 Role: Service Entry Point, User Acquisition, Monetization Gateway Tech Stack: React (Vite), Tailwind CSS, Framer Motion (Animation), Supabase Auth, Lemon Squeezy (Payments)

graph TD
    A[Landing Page (/)] -->|Click 'Login'| B(Login Modal)
    A -->|Click 'Get Started'| C(Signup Page)
    A -->|Scroll Down| D(Features & Demo)
    A -->|Scroll Down| E(Pricing Section)
    E -->|Select Plan| F(Lemon Squeezy Checkout Overlay)
    B -->|Success| G[User Dashboard]
    C -->|Success| H[Onboarding Step 1]
    F -->|Payment Success| G

    2. Landing Page UI Specification (랜딩 페이지)
2.1 Global Navigation Bar (GNB)
Position: Fixed Top (Sticky)

Height: 64px

Background: bg-white/80 (Backdrop Blur effect)

Elements:

[Left] Logo: "ADC-GenAI" (SVG) - 클릭 시 최상단 스크롤.

[Center] Menu: Features, Pricing, Resources (Golden Set).

[Right] Action Area:

[Button] Log in: Ghost Variant. 클릭 시 로그인 모달 오픈.

[Button] Start for Free: Solid Blue (#007AFF). 클릭 시 회원가입 페이지 이동.


2.2 Hero Section (최상단 메인)
Layout: 2-Column (Left: Text, Right: 3D Visual)

Left (Copywriting):

H1: "Accelerate ADC Discovery from Months to Minutes."

Sub: "AI-driven simulation for Linker-Payload optimization & Toxicity prediction."

[Input Group] Email Capture:

Input: "Enter work email..."

Button: Get Golden Set (Free)

Logic: 이메일 입력 후 버튼 클릭 시 -> CRM(DB)에 저장하고 -> "골든셋 PDF 다운로드 링크"를 이메일로 자동 발송. (리드 수집용)

Right (Visual):

Interactive 3D Component: MolStar 라이브러리를 사용해 회전하는 항체-약물 접합체(ADC) 3D 모델 렌더링.

Effect: 마우스 오버 시 특정 잔기(Residue)가 형광색으로 빛남.

2.3 Pricing Section (가격 정책 UI)
Design: 3-Card Layout.

Card 1: Free (Lead Magnet)

Title: "Researcher"

Price: $0

Features: "LIV-1 Basic Analysis", "View Golden Set Library".

[Button] Sign Up Free

Card 2: Pro (Main Product) - Highlighted Border

Title: "Developer"

Price: $499 / report (or Monthly Subscription)

Features: "Full Toxicity Prediction", "Patent Check", "Linker Optimization".

[Button] Buy Credits -> Triggers Payment Overlay

Card 3: Team (Enterprise)

Title: "Organization"

Price: Custom

Features: "API Access", "Dedicated Server", "Unlimited Seats".

[Button] Contact Sales -> Opens Typeform/Email Modal

3. Authentication Flow (인증 로직)
3.1 Login Modal Component
Trigger: GNB의 Log in 버튼 클릭.

UI Elements:

Title: "Welcome back to the Lab."

[Button] Continue with Google: Google OAuth (연구원 개인 계정).

[Button] Continue with LinkedIn: LinkedIn OAuth (비즈니스 계정 - 권장).

Divider: "Or with email"

Input: Email, Password.

[Button] Sign In (Primary)

State Handling:

isLoading: 버튼 내 스피너(Spinner) 표시, 버튼 비활성화.

isError: 인풋 하단에 빨간색 텍스트 "Invalid credentials." 표시.


// Pseudo-code for Login Button Click
const handleLogin = async (provider) => {
  setLoading(true);
  try {
    // 1. Supabase Auth 호출
    const { user, session, error } = await supabase.auth.signInWithOAuth({
      provider: 'linkedin',
    });

    if (error) throw error;

    // 2. 로그인 성공 시, 백엔드에 사용자 정보 동기화 및 크레딧 확인
    const userInfo = await api.get('/user/profile');
    
    // 3. 대시보드로 이동
    navigate('/dashboard');
  } catch (err) {
    toast.error("Login failed: " + err.message);
  } finally {
    setLoading(false);
  }
};


4. Payment Integration (결제 로직)
Provider: Lemon Squeezy (Global MoR)

4.1 Purchase Flow (구매 프로세스)
User Action: Pricing 섹션에서 [Button] Subscribe Pro 클릭.

Frontend Logic:

사용자가 로그인 상태인지 체크 (isLoggedIn).

비로그인 시 -> 3.1 Login Modal 먼저 띄움.

로그인 시 -> Lemon Squeezy의 Checkout URL을 호출하되, Overlay(팝업) 모드로 실행.

Checkout Overlay (Lemon Squeezy 제공 UI):

사용자는 여기서 카드 번호를 입력. (PG사 화면 이탈 없음)

Pre-filled Data: 로그인한 사용자의 이메일을 미리 채워줌 (&checkout[email]=user@email.com).

4.2 Post-Payment Logic (결제 완료 후 처리)
Scenario: 사용자가 결제를 완료하고 팝업을 닫음.

Frontend Polling:

결제가 서버(Webhook)에 반영될 때까지 약 3~5초 딜레이가 있을 수 있음.

결제창 닫힘 이벤트 감지 -> 로딩 스피너 화면 표시 ("Updating your license...") -> 2초마다 GET /user/subscription 호출.

응답이 plan: 'pro'로 바뀌면 -> "결제 성공! Pro 기능을 잠금 해제했습니다." (Confetti Effect 🎉) 띄우고 대시보드로 이동.


5. API Requirements (백엔드 개발자 전달용)
프론트엔드에서 인증 및 랜딩 페이지 기능을 구현하기 위해 백엔드에 필요한 API 명세입니다.

Method,Endpoint,Description,Request Body,Response Example
POST,/api/auth/lead-magnet,랜딩페이지 이메일 수집,"{ ""email"": ""..."" }","{ ""status"": ""sent"", ""msg"": ""Check inbox"" }"
GET,/api/user/profile,사용자 정보 및 크레딧 조회,(Header: Bearer Token),"{ ""id"": ""u_1"", ""credits"": 50, ""plan"": ""free"" }"
POST,/api/payment/create-checkout,레몬스퀴지 결제 링크 생성,"{ ""plan_id"": ""pro_monthly"" }","{ ""checkout_url"": ""https://lemon..."" }"


6. Development Checklist (작업 순서)
[ ] Shadcn/UI 설치 및 테마(폰트, 컬러) 설정.

[ ] Landing Page 퍼블리싱 (반응형 모바일 뷰 포함).

[ ] Supabase Auth 연동 (Google/LinkedIn).

[ ] Lemon Squeezy 테스트 모드 연동 및 결제 팝업 띄우기 확인.

[ ] Onboarding (회원가입 직후 "어떤 타겟 연구 중이세요?" 묻는 모달) 구현.