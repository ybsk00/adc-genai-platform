import { useState, useEffect, useCallback, useRef } from 'react'
import { useParams, Link } from 'react-router-dom'
import { motion } from 'framer-motion'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import {
    CheckCircle,
    Loader2,
    Clock,
    Download,
    ExternalLink,
    Dna,
    Shield,
    FileSearch,
    Users,
    Stethoscope,
    FileText,
    AlertTriangle,
    RefreshCw
} from 'lucide-react'

/**
 * [Dev Note: Adaptive Polling]
 * 초기엔 3초, 5분 뒤엔 10초로 간격을 늘려 서버 부하 감소
 */
function useAdaptivePolling(callback: () => void, isPolling: boolean) {
    const startTimeRef = useRef(Date.now())

    useEffect(() => {
        if (!isPolling) return

        const getInterval = () => {
            const elapsed = Date.now() - startTimeRef.current
            const fiveMinutes = 5 * 60 * 1000
            return elapsed > fiveMinutes ? 10000 : 3000 // 5분 후 10초, 이전엔 3초
        }

        const tick = () => {
            callback()
            timeoutId = setTimeout(tick, getInterval())
        }

        let timeoutId = setTimeout(tick, getInterval())
        return () => clearTimeout(timeoutId)
    }, [callback, isPolling])
}

// Agent 타입 정의
type AgentStatus = 'pending' | 'running' | 'done' | 'error'

interface AgentInfo {
    id: string
    name: string
    icon: React.ElementType
    description: string
    status: AgentStatus
    errorMessage?: string
}

const initialAgents: Omit<AgentInfo, 'status' | 'errorMessage'>[] = [
    { id: 'structure', name: 'Structure Agent', icon: Dna, description: '3D 구조 분석' },
    { id: 'toxicology', name: 'Toxicology Agent', icon: Shield, description: '독성 예측' },
    { id: 'patent', name: 'Patent Agent', icon: FileSearch, description: '특허 분석' },
    { id: 'competitor', name: 'Competitor Agent', icon: Users, description: '경쟁사 분석' },
    { id: 'clinical', name: 'Clinical Agent', icon: Stethoscope, description: '임상 설계' },
    { id: 'report', name: 'Report Agent', icon: FileText, description: '리포트 생성' },
]

export function ResultViewer() {
    const { jobId } = useParams()
    const [status, setStatus] = useState<'processing' | 'completed' | 'partial'>('processing')
    const [progress, setProgress] = useState(0)
    const [agents, setAgents] = useState<AgentInfo[]>(
        initialAgents.map(a => ({ ...a, status: 'pending' as const }))
    )
    const [isDownloading, setIsDownloading] = useState(false)

    // 폴링 콜백
    const pollStatus = useCallback(() => {
        // TODO: 실제 API 호출로 대체
        // const response = await fetch(`/api/jobs/${jobId}/status`)
        setProgress(prev => Math.min(prev + 5, 100))
    }, [])

    // Adaptive Polling 적용
    useAdaptivePolling(pollStatus, status === 'processing')

    // Mock: Progress에 따른 Agent 상태 업데이트
    useEffect(() => {
        if (progress >= 100) {
            // 시뮬레이션: Competitor Agent가 실패한 경우
            const hasError = Math.random() > 0.7 // 30% 확률로 부분 실패

            setAgents(prev => prev.map(agent => {
                if (agent.id === 'competitor' && hasError) {
                    return {
                        ...agent,
                        status: 'error' as const,
                        errorMessage: '네트워크 오류 - 외부 API 응답 없음'
                    }
                }
                return { ...agent, status: 'done' as const }
            }))

            setStatus(hasError ? 'partial' : 'completed')
        } else {
            // 진행 중 상태 업데이트
            const thresholds = [
                { progress: 15, agentId: 'structure' },
                { progress: 35, agentId: 'toxicology' },
                { progress: 55, agentId: 'patent' },
                { progress: 70, agentId: 'competitor' },
                { progress: 85, agentId: 'clinical' },
                { progress: 100, agentId: 'report' },
            ]

            setAgents(prev => prev.map(agent => {
                const threshold = thresholds.find(t => t.agentId === agent.id)
                if (!threshold) return agent

                if (progress >= threshold.progress) {
                    return { ...agent, status: 'done' as const }
                } else if (progress >= threshold.progress - 15) {
                    return { ...agent, status: 'running' as const }
                }
                return agent
            }))
        }
    }, [progress])

    // 부분 실패한 Agent 재시도
    const handleRetryAgent = async (agentId: string) => {
        setAgents(prev => prev.map(agent =>
            agent.id === agentId
                ? { ...agent, status: 'running' as const, errorMessage: undefined }
                : agent
        ))

        // TODO: API 호출로 특정 Agent 재실행
        await new Promise(resolve => setTimeout(resolve, 2000))

        setAgents(prev => prev.map(agent =>
            agent.id === agentId
                ? { ...agent, status: 'done' as const }
                : agent
        ))

        // 모든 Agent가 완료되면 상태 업데이트
        setAgents(prev => {
            const allDone = prev.every(a => a.status === 'done')
            if (allDone) setStatus('completed')
            return prev
        })
    }

    /**
     * [Dev Note: PDF 다운로드 보안]
     * Pre-signed URL (유효기간 5분)을 받아 다운로드
     * 링크 유출되어도 5분 후 만료됨
     */
    const handleDownloadPDF = async () => {
        setIsDownloading(true)
        try {
            // TODO: API 호출로 Pre-signed URL 획득
            // const response = await fetch(`/api/jobs/${jobId}/report-url`)
            // const { url } = await response.json()

            await new Promise(resolve => setTimeout(resolve, 1000))
            const mockPresignedUrl = `https://s3.amazonaws.com/adc-reports/${jobId}.pdf?X-Amz-Expires=300&X-Amz-Signature=xxx`

            // 다운로드 트리거
            window.open(mockPresignedUrl, '_blank')
        } finally {
            setIsDownloading(false)
        }
    }

    // 처리 중 화면
    if (status === 'processing') {
        return (
            <div className="max-w-3xl mx-auto">
                <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                >
                    <Card>
                        <CardHeader>
                            <CardTitle className="flex items-center gap-2">
                                <Loader2 className="w-5 h-5 animate-spin text-[#007AFF]" />
                                Simulation in Progress
                            </CardTitle>
                            <CardDescription>Job ID: {jobId}</CardDescription>
                        </CardHeader>
                        <CardContent>
                            {/* Progress Bar */}
                            <div className="mb-8">
                                <div className="flex items-center justify-between mb-2">
                                    <span className="text-sm text-gray-500">Overall Progress</span>
                                    <span className="text-sm font-medium text-gray-900">{progress}%</span>
                                </div>
                                <div className="h-3 bg-gray-100 rounded-full overflow-hidden">
                                    <motion.div
                                        className="h-full bg-gradient-to-r from-[#007AFF] to-purple-500"
                                        initial={{ width: 0 }}
                                        animate={{ width: `${progress}%` }}
                                        transition={{ duration: 0.3 }}
                                    />
                                </div>
                            </div>

                            {/* Agent Status */}
                            <div className="space-y-3">
                                {agents.map((agent) => (
                                    <AgentStatusRow key={agent.id} agent={agent} onRetry={handleRetryAgent} />
                                ))}
                            </div>
                        </CardContent>
                    </Card>
                </motion.div>
            </div>
        )
    }

    // 완료 또는 부분 완료 화면
    return (
        <div className="max-w-4xl mx-auto space-y-6">
            {/* 부분 실패 경고 */}
            {status === 'partial' && (
                <motion.div
                    initial={{ opacity: 0, y: -10 }}
                    animate={{ opacity: 1, y: 0 }}
                >
                    <Card className="border-amber-200 bg-amber-50">
                        <CardContent className="p-4 flex items-center gap-3">
                            <AlertTriangle className="w-5 h-5 text-amber-600" />
                            <div className="flex-1">
                                <p className="font-medium text-amber-800">일부 분석이 완료되지 않았습니다</p>
                                <p className="text-sm text-amber-600">아래에서 실패한 Agent를 재시도할 수 있습니다. 크레딧은 차감되지 않습니다.</p>
                            </div>
                        </CardContent>
                    </Card>
                </motion.div>
            )}

            {/* Summary Card */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
            >
                <Card className="bg-gradient-to-r from-blue-50 to-purple-50 border-none">
                    <CardContent className="p-8">
                        <div className="flex items-center justify-between">
                            <div>
                                <p className="text-sm text-gray-500 mb-1">Final Assessment</p>
                                <div className="flex items-center gap-4">
                                    <span className="text-6xl font-bold text-[#007AFF]">B+</span>
                                    <div>
                                        <Badge className="bg-amber-100 text-amber-700 border-amber-200 text-lg px-3 py-1">
                                            Conditional Go
                                        </Badge>
                                        <p className="text-gray-600 mt-2 max-w-md">
                                            효능은 탁월하나 독성 리스크가 발견됨.
                                            추가 in-vitro 테스트 권장.
                                        </p>
                                    </div>
                                </div>
                            </div>
                            <div className="flex flex-col gap-2">
                                <Button
                                    className="bg-[#007AFF] hover:bg-[#0056b3]"
                                    onClick={handleDownloadPDF}
                                    disabled={isDownloading}
                                >
                                    {isDownloading ? (
                                        <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                                    ) : (
                                        <Download className="w-4 h-4 mr-2" />
                                    )}
                                    Download PDF Report
                                </Button>
                                <Button variant="outline">
                                    <ExternalLink className="w-4 h-4 mr-2" />
                                    View 3D Structure
                                </Button>
                            </div>
                        </div>
                    </CardContent>
                </Card>
            </motion.div>

            {/* 부분 실패 Agent 표시 */}
            {status === 'partial' && (
                <Card>
                    <CardHeader>
                        <CardTitle>Agent Status</CardTitle>
                    </CardHeader>
                    <CardContent>
                        <div className="space-y-3">
                            {agents.map((agent) => (
                                <AgentStatusRow key={agent.id} agent={agent} onRetry={handleRetryAgent} />
                            ))}
                        </div>
                    </CardContent>
                </Card>
            )}

            {/* 3D Viewer Placeholder */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.1 }}
            >
                <Card>
                    <CardHeader>
                        <CardTitle>3D Structure View</CardTitle>
                        <CardDescription>MolStar 3D 뷰어 - 마우스로 회전, 확대 가능</CardDescription>
                    </CardHeader>
                    <CardContent>
                        <div className="aspect-video bg-gradient-to-br from-gray-900 to-gray-800 rounded-lg flex items-center justify-center">
                            <div className="text-center text-white">
                                <div className="w-24 h-24 mx-auto mb-4 rounded-full bg-gradient-to-br from-[#007AFF] to-purple-500 flex items-center justify-center">
                                    <span className="text-4xl">🧬</span>
                                </div>
                                <p className="text-gray-400">Interactive 3D ADC Model</p>
                                <p className="text-gray-500 text-sm">(MolStar Integration)</p>
                            </div>
                        </div>
                    </CardContent>
                </Card>
            </motion.div>

            {/* Radar Chart & Key Findings */}
            <div className="grid md:grid-cols-2 gap-6">
                <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    transition={{ delay: 0.2 }}
                >
                    <Card className="h-full">
                        <CardHeader>
                            <CardTitle>Analysis Scores</CardTitle>
                        </CardHeader>
                        <CardContent>
                            <div className="space-y-4">
                                {[
                                    { label: '효능 (Efficacy)', score: 85, color: 'bg-green-500' },
                                    { label: '독성 (Toxicity)', score: 60, color: 'bg-red-500' },
                                    { label: '물성 (Properties)', score: 75, color: 'bg-blue-500' },
                                    { label: '특허 (Patent)', score: 90, color: 'bg-purple-500' },
                                    { label: '사업성 (Market)', score: 70, color: 'bg-amber-500' },
                                ].map((item) => (
                                    <div key={item.label}>
                                        <div className="flex justify-between text-sm mb-1">
                                            <span className="text-gray-600">{item.label}</span>
                                            <span className="font-medium">{item.score}/100</span>
                                        </div>
                                        <div className="h-2 bg-gray-100 rounded-full overflow-hidden">
                                            <motion.div
                                                className={`h-full ${item.color}`}
                                                initial={{ width: 0 }}
                                                animate={{ width: `${item.score}%` }}
                                                transition={{ duration: 0.5, delay: 0.3 }}
                                            />
                                        </div>
                                    </div>
                                ))}
                            </div>
                        </CardContent>
                    </Card>
                </motion.div>

                <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    transition={{ delay: 0.3 }}
                >
                    <Card className="h-full">
                        <CardHeader>
                            <CardTitle>Key Findings</CardTitle>
                        </CardHeader>
                        <CardContent>
                            <div className="space-y-4">
                                <div className="flex items-start gap-3">
                                    <div className="p-1.5 bg-red-100 rounded-lg">
                                        <Shield className="w-4 h-4 text-red-600" />
                                    </div>
                                    <div>
                                        <p className="font-medium text-gray-900">독성 리스크 발견</p>
                                        <p className="text-sm text-gray-500">Neutropenia 위험 (LogP: 3.8)</p>
                                    </div>
                                </div>
                                <div className="flex items-start gap-3">
                                    <div className="p-1.5 bg-green-100 rounded-lg">
                                        <FileSearch className="w-4 h-4 text-green-600" />
                                    </div>
                                    <div>
                                        <p className="font-medium text-gray-900">특허 안전</p>
                                        <p className="text-sm text-gray-500">주요 특허 만료됨 (2024)</p>
                                    </div>
                                </div>
                                <div className="flex items-start gap-3">
                                    <div className="p-1.5 bg-blue-100 rounded-lg">
                                        <Users className="w-4 h-4 text-blue-600" />
                                    </div>
                                    <div>
                                        <p className="font-medium text-gray-900">경쟁사 현황</p>
                                        <p className="text-sm text-gray-500">Merck (Phase 2), Seagen (Phase 1)</p>
                                    </div>
                                </div>
                            </div>
                        </CardContent>
                    </Card>
                </motion.div>
            </div>

            {/* Back Button */}
            <div className="flex justify-center pt-4">
                <Button variant="outline" asChild>
                    <Link to="/dashboard">← Back to Dashboard</Link>
                </Button>
            </div>
        </div>
    )
}

// Agent 상태 표시 컴포넌트
function AgentStatusRow({
    agent,
    onRetry
}: {
    agent: AgentInfo
    onRetry: (id: string) => void
}) {
    const Icon = agent.icon

    const bgColor = {
        pending: 'bg-gray-50',
        running: 'bg-blue-50',
        done: 'bg-green-50',
        error: 'bg-red-50',
    }[agent.status]

    const iconBgColor = {
        pending: 'bg-gray-100',
        running: 'bg-blue-100',
        done: 'bg-green-100',
        error: 'bg-red-100',
    }[agent.status]

    const iconColor = {
        pending: 'text-gray-400',
        running: 'text-blue-600',
        done: 'text-green-600',
        error: 'text-red-600',
    }[agent.status]

    return (
        <div className={`flex items-center gap-4 p-3 rounded-lg transition-colors ${bgColor}`}>
            <div className={`p-2 rounded-lg ${iconBgColor}`}>
                <Icon className={`w-5 h-5 ${iconColor}`} />
            </div>
            <div className="flex-1">
                <p className="font-medium text-gray-900">{agent.name}</p>
                <p className="text-sm text-gray-500">
                    {agent.errorMessage || agent.description}
                </p>
            </div>
            <div className="flex items-center gap-2">
                {agent.status === 'done' && (
                    <Badge className="bg-green-100 text-green-700 border-green-200">
                        <CheckCircle className="w-3 h-3 mr-1" />
                        Done
                    </Badge>
                )}
                {agent.status === 'running' && (
                    <Badge className="bg-blue-100 text-blue-700 border-blue-200">
                        <Loader2 className="w-3 h-3 mr-1 animate-spin" />
                        Running
                    </Badge>
                )}
                {agent.status === 'pending' && (
                    <Badge variant="outline" className="text-gray-400">
                        <Clock className="w-3 h-3 mr-1" />
                        Pending
                    </Badge>
                )}
                {agent.status === 'error' && (
                    <>
                        <Badge className="bg-red-100 text-red-700 border-red-200">
                            <AlertTriangle className="w-3 h-3 mr-1" />
                            Error
                        </Badge>
                        <Button
                            size="sm"
                            variant="outline"
                            className="h-7 text-xs"
                            onClick={() => onRetry(agent.id)}
                        >
                            <RefreshCw className="w-3 h-3 mr-1" />
                            Retry
                        </Button>
                    </>
                )}
            </div>
        </div>
    )
}
