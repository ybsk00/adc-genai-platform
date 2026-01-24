서비스의 뼈대이자, **데이터 자산(Data Asset)**이 쌓이는 창고인 06_Database_Schema.md 문서를 작성해 드립니다.

이 문서는 Supabase (PostgreSQL) 환경을 기준으로 작성되었으며, 단순한 데이터 저장을 넘어 **RAG(벡터 검색)**와 **보안(RLS)**까지 고려한 SaaS 최적화 스키마입니다.

06. Database Schema & Vector Store Design
Document ID: DB-01 Role: Data Storage, Vector Search Engine, User Management Tech Stack: Supabase (PostgreSQL 15+), pgvector (Vector Extension)

1. Entity-Relationship Diagram (ERD)
데이터 간의 관계를 한눈에 파악하는 지도입니다.

코드 스니펫

erDiagram
    organizations ||--|{ profiles : "has members"
    profiles ||--|{ projects : "creates"
    profiles ||--|{ transactions : "makes"
    
    projects {
        uuid id PK
        uuid user_id FK
        jsonb input_data
        jsonb result_data
        string status
    }
    
    golden_set_library ||--o{ golden_set_embeddings : "has vectors"
    
    golden_set_library {
        uuid id PK
        string name
        string category
        jsonb properties
        string patent_status
    }
    
    agent_prompts {
        string agent_id PK
        text content
        int version
    }

    2. Core Tables Specification (핵심 테이블 명세)
2.1 User & Organization (회원 및 조직)
Supabase의 기본 인증(auth.users)과 연동되는 프로필 테이블입니다.

👥 profiles (Public User Data)
Description: 사용자 상세 정보 및 크레딧 잔액.

Primary Key: id (references auth.users.id)

Column Name,Data Type,Constraint,Description
id,uuid,"PK, FK",Supabase Auth ID와 1:1 매칭
email,text,Unique,이메일 주소
full_name,text,,사용자 실명
org_id,uuid,FK,소속 조직 ID (Optional)
role,text,,"'user', 'admin', 'viewer'"
credits,int,Default 0,현재 보유 크레딧 (화폐)
created_at,timestamptz,,가입일

Column Name,Data Type,Description
id,uuid,PK
name,text,"조직명 (e.g., ""Syngene Team A"")"
domain,text,"자동 가입 도메인 (e.g., @syngene.com)"
subscription_tier,text,"'free', 'pro', 'enterprise'"


2.2 Simulation Engine (프로젝트 및 결과)
사용자가 돌린 시뮬레이션의 **입력값(레시피)**과 **결과값(요리)**을 저장합니다.

🧪 projects (Simulation Jobs)
Description: 하나의 시뮬레이션 작업 단위.



Column Name,Data Type,Description
id,uuid,PK (Job ID)
user_id,uuid,FK (작성자)
title,text,"프로젝트명 (e.g., ""LIV-1 Test 01"")"
status,text,"'queued', 'processing', 'completed', 'failed'"
input_data,jsonb,"입력 데이터 (seq, smiles, mode)"
result_data,jsonb,AI 분석 결과 (6-Agent Output JSON)
report_url,text,최종 PDF 다운로드 링크 (S3 URL)
error_log,text,실패 시 에러 메시지 저장
created_at,timestamptz,생성 시간


JSONB를 쓰는 이유: 바이오 데이터(구조 정보, 독성 수치 등)는 항목이 너무 다양해서 컬럼을 고정하면 나중에 수정하기 힘듭니다. jsonb로 유연하게 저장합니다.

2.3 Knowledge Base (RAG & Vector)
우리 서비스의 핵심 자산인 **Golden Set(족보)**과 이를 AI가 검색하기 위한 벡터(Vector) 데이터입니다.


Column Name,Data Type,Description
id,uuid,PK
name,text,"물질 이름 (e.g., ""MMAE"", ""Trastuzumab"")"
category,text,"'antibody', 'linker', 'payload', 'drug'"
structure_code,text,FASTA Sequence 또는 SMILES Code
properties,jsonb,"물성 정보 (logP, MW, clinical_stage 등)"
description,text,텍스트 설명 (RAG 임베딩용 원문)
source,text,"출처 (e.g., ""ClinicalTrials.gov"", ""Patent US-123"")"


golden_set_embeddings (Vector Store)
Description: golden_set_library의 텍스트 설명을 **숫자(Vector)**로 변환한 테이블.

  similarity float
)
language plpgsql
as $$
begin
  return query
  select
    golden_set_embeddings.source_id,
    golden_set_embeddings.chunk_content,
    1 - (golden_set_embeddings.embedding <=> query_embedding) as similarity
  from golden_set_embeddings
  where 1 - (golden_set_embeddings.embedding <=> query_embedding) > match_threshold
  order by similarity desc
  limit match_count;
end;
$$;


💡 작성자의 코멘트
이 스키마 설계의 핵심은 **"유연함(Flexibility)"**입니다.

projects.result_data (JSONB): AI 에이전트가 4명에서 6명, 10명으로 늘어나도 DB 스키마를 뜯어고칠 필요가 없습니다. 그냥 JSON에 필드만 추가하면 됩니다.

golden_set_embeddings: RAG 시스템을 위한 준비가 DB 레벨에서 완벽하게 끝났습니다.