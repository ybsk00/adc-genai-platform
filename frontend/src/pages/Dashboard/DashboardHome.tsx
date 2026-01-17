import { useState, useEffect } from 'react'
import { motion } from 'framer-motion'
import { Link } from 'react-router-dom'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import { Skeleton } from '@/components/ui/skeleton'
import {
    Plus,
    Library,
    Clock,
    CheckCircle,
    XCircle,
    Loader2
} from 'lucide-react'
import { RealTimeFeed } from '@/components/dashboard/RealTimeFeed'

// Mock data for recent simulations
const recentSimulations = [
    { id: 'job_1', name: 'LIV-1 + MMAE Test', status: 'completed', grade: 'B+', date: '2026-01-17' },
    { id: 'job_2', name: 'HER2 + DXd Analysis', status: 'processing', progress: 65, date: '2026-01-17' },
    { id: 'job_3', name: 'TROP2 + SN-38', status: 'failed', error: 'Timeout', date: '2026-01-16' },
    { id: 'job_4', name: 'EGFR + MMAF', status: 'completed', grade: 'A', date: '2026-01-15' },
]

const statusConfig = {
    completed: { icon: CheckCircle, color: 'text-green-500', bg: 'bg-green-50' },
    processing: { icon: Loader2, color: 'text-blue-500', bg: 'bg-blue-50' },
    failed: { icon: XCircle, color: 'text-red-500', bg: 'bg-red-50' },
}

/**
 * [Dev Note: Skeleton Loading]
 * 스피너 대신 스켈레톤 UI로 체감 로딩 속도 향상
 */
function DashboardSkeleton() {
    return (
        <div className="space-y-6">
            {/* Welcome Skeleton */}
            <div>
                <Skeleton className="h-8 w-80 mb-2" />
                <Skeleton className="h-5 w-60" />
            </div>

            {/* Quick Actions Skeleton */}
            <div className="grid md:grid-cols-2 gap-4">
                <Skeleton className="h-24 rounded-xl" />
                <Skeleton className="h-24 rounded-xl" />
            </div>

            {/* Recent Activity Skeleton */}
            <Card>
                <CardHeader>
                    <Skeleton className="h-6 w-40" />
                    <Skeleton className="h-4 w-32" />
                </CardHeader>
                <CardContent>
                    <div className="space-y-3">
                        {[1, 2, 3, 4].map((i) => (
                            <div key={i} className="flex items-center gap-3 p-3">
                                <Skeleton className="w-10 h-10 rounded-lg" />
                                <div className="flex-1">
                                    <Skeleton className="h-5 w-40 mb-1" />
                                    <Skeleton className="h-4 w-24" />
                                </div>
                                <Skeleton className="h-6 w-20 rounded-full" />
                            </div>
                        ))}
                    </div>
                </CardContent>
            </Card>

            {/* Trend Feed Skeleton */}
            <Card>
                <CardHeader>
                    <Skeleton className="h-6 w-40" />
                    <Skeleton className="h-4 w-48" />
                </CardHeader>
                <CardContent>
                    <div className="space-y-4">
                        <Skeleton className="h-20 rounded-lg" />
                        <Skeleton className="h-20 rounded-lg" />
                    </div>
                </CardContent>
            </Card>
        </div>
    )
}

export function DashboardHome() {
    const [isLoading, setIsLoading] = useState(true)
    const [userName] = useState('Dr. Kim')

    // 데이터 로딩 시뮬레이션
    useEffect(() => {
        // TODO: 실제 API 호출로 대체
        const timer = setTimeout(() => {
            setIsLoading(false)
        }, 1500)
        return () => clearTimeout(timer)
    }, [])

    // 스켈레톤 로딩 표시
    if (isLoading) {
        return <DashboardSkeleton />
    }

    return (
        <div className="space-y-6">
            {/* Welcome Message */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
            >
                <h1 className="text-2xl font-bold text-gray-900">
                    Good Morning, {userName}. Ready to discover?
                </h1>
                <p className="text-gray-500 mt-1">
                    오늘도 새로운 ADC 후보를 분석해보세요.
                </p>
            </motion.div>

            {/* Quick Actions */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.1 }}
                className="grid md:grid-cols-2 gap-4"
            >
                <Card className="border-2 border-dashed border-[#007AFF] bg-blue-50/50 hover:bg-blue-50 transition-colors cursor-pointer">
                    <Link to="/dashboard/builder">
                        <CardContent className="flex items-center gap-4 p-6">
                            <div className="w-14 h-14 rounded-xl bg-[#007AFF] flex items-center justify-center">
                                <Plus className="w-7 h-7 text-white" />
                            </div>
                            <div>
                                <h3 className="font-semibold text-gray-900 text-lg">New Simulation</h3>
                                <p className="text-gray-500 text-sm">AI 기반 ADC 분석 시작하기</p>
                            </div>
                        </CardContent>
                    </Link>
                </Card>

                <Card className="hover:shadow-md transition-shadow cursor-pointer">
                    <Link to="/dashboard/library">
                        <CardContent className="flex items-center gap-4 p-6">
                            <div className="w-14 h-14 rounded-xl bg-gray-100 flex items-center justify-center">
                                <Library className="w-7 h-7 text-gray-600" />
                            </div>
                            <div>
                                <h3 className="font-semibold text-gray-900 text-lg">Browse Golden Set</h3>
                                <p className="text-gray-500 text-sm">FDA 승인 ADC 데이터베이스 검색</p>
                            </div>
                        </CardContent>
                    </Link>
                </Card>
            </motion.div>

            {/* Recent Activity */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.2 }}
            >
                <Card>
                    <CardHeader>
                        <CardTitle className="flex items-center gap-2">
                            <Clock className="w-5 h-5" />
                            Recent Activity
                        </CardTitle>
                        <CardDescription>최근 시뮬레이션 현황</CardDescription>
                    </CardHeader>
                    <CardContent>
                        <div className="space-y-3">
                            {recentSimulations.map((sim) => {
                                const config = statusConfig[sim.status as keyof typeof statusConfig]
                                const StatusIcon = config.icon
                                return (
                                    <Link
                                        key={sim.id}
                                        to={`/dashboard/result/${sim.id}`}
                                        className="flex items-center justify-between p-3 rounded-lg hover:bg-gray-50 transition-colors"
                                    >
                                        <div className="flex items-center gap-3">
                                            <div className={`p-2 rounded-lg ${config.bg}`}>
                                                <StatusIcon className={`w-4 h-4 ${config.color} ${sim.status === 'processing' ? 'animate-spin' : ''}`} />
                                            </div>
                                            <div>
                                                <p className="font-medium text-gray-900">{sim.name}</p>
                                                <p className="text-sm text-gray-500">{sim.date}</p>
                                            </div>
                                        </div>
                                        <div>
                                            {sim.status === 'completed' && (
                                                <Badge variant="outline" className="border-green-200 text-green-700 bg-green-50">
                                                    Grade: {sim.grade}
                                                </Badge>
                                            )}
                                            {sim.status === 'processing' && (
                                                <Badge variant="outline" className="border-blue-200 text-blue-700 bg-blue-50">
                                                    {sim.progress}%
                                                </Badge>
                                            )}
                                            {sim.status === 'failed' && (
                                                <Badge variant="outline" className="border-red-200 text-red-700 bg-red-50">
                                                    {sim.error}
                                                </Badge>
                                            )}
                                        </div>
                                    </Link>
                                )
                            })}
                        </div>
                    </CardContent>
                </Card>
            </motion.div>

            {/* ADC Trend Feed */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.3 }}
            >
                <RealTimeFeed />
            </motion.div>
        </div>
    )
}
