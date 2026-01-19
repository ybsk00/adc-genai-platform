/**
 * API Client - 백엔드 통신을 위한 클라이언트
 * 
 * 환경변수 VITE_API_BASE_URL을 사용하여 백엔드 URL을 설정합니다.
 * 개발: http://localhost:8000
 * 프로덕션: https://adc-backend-962229188169.asia-northeast3.run.app
 */

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'https://adc-backend-962229188169.asia-northeast3.run.app'


interface ApiResponse<T> {
    data: T | null
    error: string | null
}

/**
 * API 요청을 수행하는 기본 함수
 */
async function apiRequest<T>(
    endpoint: string,
    options: RequestInit = {}
): Promise<ApiResponse<T>> {
    const url = `${API_BASE_URL}${endpoint}`

    try {
        const response = await fetch(url, {
            ...options,
            headers: {
                'Content-Type': 'application/json',
                ...options.headers,
            },
        })

        if (!response.ok) {
            const error = await response.json().catch(() => ({ detail: response.statusText }))
            return { data: null, error: error.detail || `Error: ${response.status}` }
        }

        const data = await response.json()
        return { data, error: null }
    } catch (error) {
        return { data: null, error: error instanceof Error ? error.message : 'Network error' }
    }
}

/**
 * 인증된 API 요청 (Supabase 토큰 포함)
 */
async function apiRequestWithAuth<T>(
    endpoint: string,
    token: string,
    options: RequestInit = {}
): Promise<ApiResponse<T>> {
    return apiRequest<T>(endpoint, {
        ...options,
        headers: {
            ...options.headers,
            'Authorization': `Bearer ${token}`,
        },
    })
}

// === API Endpoints ===

// Health Check
export const healthCheck = () => apiRequest('/health')

// Jobs API
export const createJob = (data: object, token: string) =>
    apiRequestWithAuth('/api/jobs', token, { method: 'POST', body: JSON.stringify(data) })

export const getJob = (jobId: string, token: string) =>
    apiRequestWithAuth(`/api/jobs/${jobId}`, token)

export const getJobStatus = (jobId: string, token: string) =>
    apiRequestWithAuth(`/api/jobs/${jobId}/status`, token)

export const getJobReport = (jobId: string, token: string) =>
    apiRequestWithAuth(`/api/jobs/${jobId}/report-url`, token)

// Library API
export const getGoldenSet = () => apiRequest('/api/library/goldenset')

export const searchGoldenSet = (query: string) =>
    apiRequest(`/api/library/search?q=${encodeURIComponent(query)}`)

// Admin API
export const adminGrantCredits = (userId: string, amount: number, token: string) =>
    apiRequestWithAuth('/api/admin/grant-credits', token, {
        method: 'POST',
        body: JSON.stringify({ user_id: userId, amount }),
    })

// RAG Upload API
export const uploadDocument = (file: File, metadata: object, token: string) => {
    const formData = new FormData()
    formData.append('file', file)
    formData.append('metadata', JSON.stringify(metadata))

    return apiRequestWithAuth<{ id: string; status: string }>('/api/admin/upload', token, {
        method: 'POST',
        body: formData,
        headers: {}, // Let browser set Content-Type for FormData
    })
}

export { API_BASE_URL, apiRequest, apiRequestWithAuth }
