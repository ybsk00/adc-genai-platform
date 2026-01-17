import { useState, useEffect, useCallback } from 'react'
import { motion } from 'framer-motion'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import {
    DollarSign,
    Users,
    Zap,
    AlertTriangle,
    TrendingUp,
    TrendingDown,
    Activity,
    CheckCircle,
    XCircle,
    Loader2
} from 'lucide-react'

// Mock KPI Data
const kpiData = [
    {
        label: 'MRR',
        value: '$12,450',
        change: '+12.5%',
        trend: 'up',
        icon: DollarSign,
        color: 'from-green-500 to-emerald-500'
    },
    {
        label: 'Total Users',
        value: '1,204',
        change: '+8.2%',
        trend: 'up',
        icon: Users,
        color: 'from-blue-500 to-cyan-500'
    },
    {
        label: 'Active Simulations',
        value: '45',
        change: '+23.1%',
        trend: 'up',
        icon: Zap,
        color: 'from-purple-500 to-pink-500'
    },
    {
        label: 'Error Rate',
        value: '0.5%',
        change: '-0.3%',
        trend: 'down',
        icon: AlertTriangle,
        color: 'from-amber-500 to-orange-500'
    },
]

// Mock Recent Activity
const recentActivity = [
    { user: 'Dr. Kim', action: 'Started Deep Analysis', time: '2 minutes ago', type: 'simulation' },
    { user: 'Prof. Lee', action: 'Downloaded Report', time: '15 minutes ago', type: 'report' },
    { user: 'System', action: 'RAG Index Updated', time: '1 hour ago', type: 'system' },
    { user: 'Dr. Park', action: 'Upgraded to Pro', time: '3 hours ago', type: 'payment' },
    { user: 'System', action: 'Clinical Data Synced', time: '6 hours ago', type: 'system' },
]

// System services to monitor
const systemServices = [
    { id: 'api', name: 'API Server', endpoint: '/api/health' },
    { id: 'agents', name: 'AI Agents', endpoint: '/api/agents/health' },
    { id: 'rag', name: 'RAG Pipeline', endpoint: '/api/rag/health' },
    { id: 'db', name: 'Supabase DB', endpoint: '/api/db/health' },
]

type ServiceStatus = 'operational' | 'degraded' | 'down' | 'checking'

/**
 * [Dev Note: Live Health Check]
 * 1분마다 각 서비스의 /health 엔드포인트 호출
 * 200 OK → operational, 그 외 → down
 */
function useHealthCheck() {
    const [statuses, setStatuses] = useState<Record<string, ServiceStatus>>(
        Object.fromEntries(systemServices.map(s => [s.id, 'checking']))
    )

    const checkHealth = useCallback(async () => {
        for (const service of systemServices) {
            try {
                // TODO: 실제 API 호출로 대체
                // const response = await fetch(service.endpoint)
                // const isHealthy = response.ok

                // Mock: 랜덤하게 상태 시뮬레이션
                await new Promise(r => setTimeout(r, 100))
                const isHealthy = Math.random() > 0.05 // 95% 확률로 정상

                setStatuses(prev => ({
                    ...prev,
                    [service.id]: isHealthy ? 'operational' : 'down'
                }))
            } catch {
                setStatuses(prev => ({
                    ...prev,
                    [service.id]: 'down'
                }))
            }
        }
    }, [])

    useEffect(() => {
        // 초기 체크
        checkHealth()

        // 1분마다 폴링
        const interval = setInterval(checkHealth, 60000)
        return () => clearInterval(interval)
    }, [checkHealth])

    return statuses
}

export function AdminOverview() {
    const serviceStatuses = useHealthCheck()

    return (
        <div className="space-y-6">
            {/* Page Header */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
            >
                <h1 className="text-2xl font-bold text-white">Admin Overview</h1>
                <p className="text-slate-400 mt-1">실시간 플랫폼 현황 모니터링</p>
            </motion.div>

            {/* KPI Cards */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.1 }}
                className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4"
            >
                {kpiData.map((kpi) => (
                    <Card key={kpi.label} className="bg-slate-900 border-slate-800">
                        <CardContent className="p-6">
                            <div className="flex items-center justify-between">
                                <div>
                                    <p className="text-sm text-slate-400">{kpi.label}</p>
                                    <p className="text-3xl font-bold text-white mt-1">{kpi.value}</p>
                                    <div className="flex items-center gap-1 mt-2">
                                        {kpi.trend === 'up' ? (
                                            <TrendingUp className="w-4 h-4 text-green-500" />
                                        ) : (
                                            <TrendingDown className="w-4 h-4 text-red-500" />
                                        )}
                                        <span className={kpi.trend === 'up' ? 'text-green-500 text-sm' : 'text-red-500 text-sm'}>
                                            {kpi.change}
                                        </span>
                                        <span className="text-slate-500 text-sm">vs last month</span>
                                    </div>
                                </div>
                                <div className={`w-12 h-12 rounded-xl bg-gradient-to-br ${kpi.color} flex items-center justify-center`}>
                                    <kpi.icon className="w-6 h-6 text-white" />
                                </div>
                            </div>
                        </CardContent>
                    </Card>
                ))}
            </motion.div>

            {/* Charts Row */}
            <div className="grid lg:grid-cols-2 gap-6">
                {/* Usage Chart Placeholder */}
                <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    transition={{ delay: 0.2 }}
                >
                    <Card className="bg-slate-900 border-slate-800">
                        <CardHeader>
                            <CardTitle className="text-white flex items-center gap-2">
                                <Activity className="w-5 h-5" />
                                Daily Usage
                            </CardTitle>
                            <CardDescription className="text-slate-400">
                                지난 7일간 시뮬레이션 요청 수
                            </CardDescription>
                        </CardHeader>
                        <CardContent>
                            {/* Chart Placeholder */}
                            <div className="h-64 bg-slate-800/50 rounded-lg flex items-center justify-center border border-slate-700">
                                <div className="text-center">
                                    <div className="flex items-end justify-center gap-2 mb-4">
                                        {[40, 65, 45, 80, 55, 90, 70].map((h, i) => (
                                            <motion.div
                                                key={i}
                                                initial={{ height: 0 }}
                                                animate={{ height: h }}
                                                transition={{ delay: 0.3 + i * 0.1 }}
                                                className="w-8 bg-gradient-to-t from-purple-500 to-pink-500 rounded-t"
                                                style={{ height: `${h}%` }}
                                            />
                                        ))}
                                    </div>
                                    <p className="text-slate-500 text-sm">Recharts 통합 예정</p>
                                </div>
                            </div>
                        </CardContent>
                    </Card>
                </motion.div>

                {/* Recent Activity */}
                <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    transition={{ delay: 0.3 }}
                >
                    <Card className="bg-slate-900 border-slate-800">
                        <CardHeader>
                            <CardTitle className="text-white">Recent Activity</CardTitle>
                            <CardDescription className="text-slate-400">
                                실시간 사용자 활동 로그
                            </CardDescription>
                        </CardHeader>
                        <CardContent>
                            <div className="space-y-4">
                                {recentActivity.map((activity, index) => (
                                    <div key={index} className="flex items-center gap-3">
                                        <div className={`w-2 h-2 rounded-full ${activity.type === 'simulation' ? 'bg-blue-500' :
                                                activity.type === 'report' ? 'bg-green-500' :
                                                    activity.type === 'payment' ? 'bg-purple-500' :
                                                        'bg-slate-500'
                                            }`} />
                                        <div className="flex-1">
                                            <p className="text-sm text-white">
                                                <span className="font-medium">{activity.user}</span>
                                                <span className="text-slate-400"> {activity.action}</span>
                                            </p>
                                            <p className="text-xs text-slate-500">{activity.time}</p>
                                        </div>
                                        <Badge
                                            variant="outline"
                                            className={`text-xs ${activity.type === 'simulation' ? 'border-blue-500/30 text-blue-400' :
                                                    activity.type === 'report' ? 'border-green-500/30 text-green-400' :
                                                        activity.type === 'payment' ? 'border-purple-500/30 text-purple-400' :
                                                            'border-slate-600 text-slate-400'
                                                }`}
                                        >
                                            {activity.type}
                                        </Badge>
                                    </div>
                                ))}
                            </div>
                        </CardContent>
                    </Card>
                </motion.div>
            </div>

            {/* System Status - Live Health Check */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.4 }}
            >
                <Card className="bg-slate-900 border-slate-800">
                    <CardHeader>
                        <div className="flex items-center justify-between">
                            <div>
                                <CardTitle className="text-white">System Status</CardTitle>
                                <CardDescription className="text-slate-400">
                                    주요 서비스 상태 (1분마다 자동 갱신)
                                </CardDescription>
                            </div>
                            <Badge variant="outline" className="border-slate-600 text-slate-400">
                                Live
                            </Badge>
                        </div>
                    </CardHeader>
                    <CardContent>
                        <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                            {systemServices.map((service) => {
                                const status = serviceStatuses[service.id]
                                return (
                                    <div
                                        key={service.id}
                                        className="flex items-center gap-3 p-3 bg-slate-800/50 rounded-lg border border-slate-700"
                                    >
                                        {status === 'checking' ? (
                                            <Loader2 className="w-4 h-4 text-slate-400 animate-spin" />
                                        ) : status === 'operational' ? (
                                            <CheckCircle className="w-4 h-4 text-green-500" />
                                        ) : (
                                            <XCircle className="w-4 h-4 text-red-500" />
                                        )}
                                        <div>
                                            <p className="text-sm font-medium text-white">{service.name}</p>
                                            <p className={`text-xs capitalize ${status === 'operational' ? 'text-green-400' :
                                                    status === 'checking' ? 'text-slate-400' :
                                                        'text-red-400'
                                                }`}>
                                                {status}
                                            </p>
                                        </div>
                                    </div>
                                )
                            })}
                        </div>
                    </CardContent>
                </Card>
            </motion.div>
        </div>
    )
}
