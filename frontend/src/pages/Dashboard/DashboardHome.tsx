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
    completed: { icon: CheckCircle, color: 'text-green-400', bg: 'bg-green-500/10' },
    processing: { icon: Loader2, color: 'text-blue-400', bg: 'bg-blue-500/10' },
    failed: { icon: XCircle, color: 'text-red-400', bg: 'bg-red-500/10' },
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
                <Skeleton className="h-8 w-80 mb-2 bg-slate-800" />
                <Skeleton className="h-5 w-60 bg-slate-800" />
            </div>

            {/* Quick Actions Skeleton */}
            <div className="grid md:grid-cols-2 gap-4">
                <Skeleton className="h-24 rounded-xl bg-slate-800" />
                <Skeleton className="h-24 rounded-xl bg-slate-800" />
            </div>

            {/* Recent Activity Skeleton */}
            <Card className="bg-slate-900 border-slate-800">
                <CardHeader>
                    <Skeleton className="h-6 w-40 bg-slate-800" />
                    <Skeleton className="h-4 w-32 bg-slate-800" />
                </CardHeader>
                <CardContent>
                    <div className="space-y-3">
                        {[1, 2, 3, 4].map((i) => (
                            <div key={i} className="flex items-center gap-3 p-3">
                                <Skeleton className="w-10 h-10 rounded-lg bg-slate-800" />
                                <div className="flex-1">
                                    <Skeleton className="h-5 w-40 mb-1 bg-slate-800" />
                                    <Skeleton className="h-4 w-24 bg-slate-800" />
                                </div>
                                <Skeleton className="h-6 w-20 rounded-full bg-slate-800" />
                            </div>
                        ))}
                    </div>
                </CardContent>
            </Card>

            {/* Trend Feed Skeleton */}
            <Card className="bg-slate-900 border-slate-800">
                <CardHeader>
                    <Skeleton className="h-6 w-40 bg-slate-800" />
                    <Skeleton className="h-4 w-48 bg-slate-800" />
                </CardHeader>
                <CardContent>
                    <div className="space-y-4">
                        <Skeleton className="h-20 rounded-lg bg-slate-800" />
                        <Skeleton className="h-20 rounded-lg bg-slate-800" />
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
                <h1 className="text-2xl font-bold text-white">
                    Good Morning, {userName}. Ready to discover?
                </h1>
                <p className="text-slate-400 mt-1">
                    Analyze new ADC candidates today.
                </p>
            </motion.div>

            {/* Quick Actions */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.1 }}
                className="grid md:grid-cols-2 gap-4"
            >
                <Card className="border-2 border-dashed border-blue-500/30 bg-blue-500/5 hover:bg-blue-500/10 transition-colors cursor-pointer">
                    <Link to="/dashboard/builder">
                        <CardContent className="flex items-center gap-4 p-6">
                            <div className="w-14 h-14 rounded-xl bg-blue-500/20 flex items-center justify-center">
                                <Plus className="w-7 h-7 text-blue-400" />
                            </div>
                            <div>
                                <h3 className="font-semibold text-white text-lg">New Simulation</h3>
                                <p className="text-slate-400 text-sm">Start AI-based ADC analysis</p>
                            </div>
                        </CardContent>
                    </Link>
                </Card>

                <Card className="bg-slate-900 border-slate-800 hover:border-slate-700 transition-colors cursor-pointer">
                    <Link to="/dashboard/library">
                        <CardContent className="flex items-center gap-4 p-6">
                            <div className="w-14 h-14 rounded-xl bg-slate-800 flex items-center justify-center">
                                <Library className="w-7 h-7 text-slate-400" />
                            </div>
                            <div>
                                <h3 className="font-semibold text-white text-lg">Browse Golden Set</h3>
                                <p className="text-slate-400 text-sm">Search FDA-approved ADC database</p>
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
                <Card className="bg-slate-900 border-slate-800">
                    <CardHeader>
                        <CardTitle className="flex items-center gap-2 text-white">
                            <Clock className="w-5 h-5 text-slate-400" />
                            Recent Activity
                        </CardTitle>
                        <CardDescription className="text-slate-400">Recent Simulation Status</CardDescription>
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
                                        className="flex items-center justify-between p-3 rounded-lg hover:bg-slate-800 transition-colors"
                                    >
                                        <div className="flex items-center gap-3">
                                            <div className={`p-2 rounded-lg ${config.bg}`}>
                                                <StatusIcon className={`w-4 h-4 ${config.color} ${sim.status === 'processing' ? 'animate-spin' : ''}`} />
                                            </div>
                                            <div>
                                                <p className="font-medium text-white">{sim.name}</p>
                                                <p className="text-sm text-slate-400">{sim.date}</p>
                                            </div>
                                        </div>
                                        <div>
                                            {sim.status === 'completed' && (
                                                <Badge variant="outline" className="border-green-500/30 text-green-400 bg-green-500/10">
                                                    Grade: {sim.grade}
                                                </Badge>
                                            )}
                                            {sim.status === 'processing' && (
                                                <Badge variant="outline" className="border-blue-500/30 text-blue-400 bg-blue-500/10">
                                                    {sim.progress}%
                                                </Badge>
                                            )}
                                            {sim.status === 'failed' && (
                                                <Badge variant="outline" className="border-red-500/30 text-red-400 bg-red-500/10">
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
