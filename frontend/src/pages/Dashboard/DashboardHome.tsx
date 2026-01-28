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
    Beaker,
    Sparkles,
    Shield,
    Factory,
    ArrowRight,
    Inbox
} from 'lucide-react'
import { RealTimeFeed } from '@/components/dashboard/RealTimeFeed'

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

            {/* Design Engine - 4대 메뉴 */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.15 }}
            >
                <div className="flex items-center justify-between mb-4">
                    <h2 className="text-lg font-semibold text-white flex items-center gap-2">
                        <Beaker className="w-5 h-5 text-purple-400" />
                        AI Design Engine
                    </h2>
                    <Badge variant="outline" className="border-purple-500/30 text-purple-400 bg-purple-500/10">
                        Multi-Agent System
                    </Badge>
                </div>
                <div className="grid md:grid-cols-4 gap-3">
                    {/* De novo Design */}
                    <Link to="/dashboard/denovo">
                        <Card className="bg-slate-900 border-slate-800 hover:border-purple-500/50 transition-all cursor-pointer group h-full">
                            <CardContent className="p-4">
                                <div className="w-10 h-10 rounded-lg bg-purple-500/20 flex items-center justify-center mb-3 group-hover:scale-110 transition-transform">
                                    <Beaker className="w-5 h-5 text-purple-400" />
                                </div>
                                <h3 className="font-medium text-white mb-1">De novo Design</h3>
                                <p className="text-xs text-slate-400 mb-3">AI-powered new ADC candidate generation</p>
                                <div className="flex items-center text-purple-400 text-xs font-medium opacity-0 group-hover:opacity-100 transition-opacity">
                                    Start Design <ArrowRight className="w-3 h-3 ml-1" />
                                </div>
                            </CardContent>
                        </Card>
                    </Link>

                    {/* Lead Optimization */}
                    <Link to="/dashboard/optimization">
                        <Card className="bg-slate-900 border-slate-800 hover:border-amber-500/50 transition-all cursor-pointer group h-full">
                            <CardContent className="p-4">
                                <div className="w-10 h-10 rounded-lg bg-amber-500/20 flex items-center justify-center mb-3 group-hover:scale-110 transition-transform">
                                    <Sparkles className="w-5 h-5 text-amber-400" />
                                </div>
                                <h3 className="font-medium text-white mb-1">Lead Optimization</h3>
                                <p className="text-xs text-slate-400 mb-3">Optimize existing compound properties</p>
                                <div className="flex items-center text-amber-400 text-xs font-medium opacity-0 group-hover:opacity-100 transition-opacity">
                                    Start Optimization <ArrowRight className="w-3 h-3 ml-1" />
                                </div>
                            </CardContent>
                        </Card>
                    </Link>

                    {/* Pre-clinical Audit */}
                    <Link to="/dashboard/audit">
                        <Card className="bg-slate-900 border-slate-800 hover:border-blue-500/50 transition-all cursor-pointer group h-full">
                            <CardContent className="p-4">
                                <div className="w-10 h-10 rounded-lg bg-blue-500/20 flex items-center justify-center mb-3 group-hover:scale-110 transition-transform">
                                    <Shield className="w-5 h-5 text-blue-400" />
                                </div>
                                <h3 className="font-medium text-white mb-1">Pre-clinical Audit</h3>
                                <p className="text-xs text-slate-400 mb-3">Safety & efficacy assessment</p>
                                <div className="flex items-center text-blue-400 text-xs font-medium opacity-0 group-hover:opacity-100 transition-opacity">
                                    Start Audit <ArrowRight className="w-3 h-3 ml-1" />
                                </div>
                            </CardContent>
                        </Card>
                    </Link>

                    {/* CMC & Sourcing */}
                    <Link to="/dashboard/cmc">
                        <Card className="bg-slate-900 border-slate-800 hover:border-emerald-500/50 transition-all cursor-pointer group h-full">
                            <CardContent className="p-4">
                                <div className="w-10 h-10 rounded-lg bg-emerald-500/20 flex items-center justify-center mb-3 group-hover:scale-110 transition-transform">
                                    <Factory className="w-5 h-5 text-emerald-400" />
                                </div>
                                <h3 className="font-medium text-white mb-1">CMC & Sourcing</h3>
                                <p className="text-xs text-slate-400 mb-3">Manufacturing & supplier analysis</p>
                                <div className="flex items-center text-emerald-400 text-xs font-medium opacity-0 group-hover:opacity-100 transition-opacity">
                                    Start Analysis <ArrowRight className="w-3 h-3 ml-1" />
                                </div>
                            </CardContent>
                        </Card>
                    </Link>
                </div>
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
                        {/* Empty State */}
                        <div className="flex flex-col items-center justify-center py-12 text-center">
                            <div className="w-16 h-16 rounded-full bg-slate-800 flex items-center justify-center mb-4">
                                <Inbox className="w-8 h-8 text-slate-500" />
                            </div>
                            <h3 className="text-lg font-medium text-white mb-2">No recent activity</h3>
                            <p className="text-slate-400 text-sm mb-4 max-w-sm">
                                Start a new ADC simulation to see your activity here.
                            </p>
                            <Link to="/dashboard/builder">
                                <Badge variant="outline" className="border-blue-500/30 text-blue-400 bg-blue-500/10 cursor-pointer hover:bg-blue-500/20 transition-colors">
                                    <Plus className="w-3 h-3 mr-1" />
                                    New Simulation
                                </Badge>
                            </Link>
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
