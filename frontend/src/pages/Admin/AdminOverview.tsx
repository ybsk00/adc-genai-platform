import { useState, useEffect, useCallback } from 'react'
import { motion } from 'framer-motion'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import {
    FileText,
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
import {
    LineChart,
    Line,
    XAxis,
    YAxis,
    CartesianGrid,
    Tooltip,
    ResponsiveContainer,
    PieChart,
    Pie,
    Cell
} from 'recharts'
import { API_BASE_URL } from '@/lib/api'

// Mock KPI Data (To be replaced with API data)
const kpiData = [
    {
        label: 'Docs Crawled Today',
        value: '1,240',
        change: '+150',
        trend: 'up',
        icon: FileText,
        color: 'from-blue-500 to-cyan-500'
    },
    {
        label: 'Simulations Ran Today',
        value: '45',
        change: '+12',
        trend: 'up',
        icon: Zap,
        color: 'from-purple-500 to-pink-500'
    },
    {
        label: 'Error Rate',
        value: '0.5%',
        change: '-0.2%',
        trend: 'down',
        icon: AlertTriangle,
        color: 'from-amber-500 to-orange-500'
    },
    {
        label: 'Total Users',
        value: '1,204',
        change: '+15',
        trend: 'up',
        icon: Users,
        color: 'from-green-500 to-emerald-500'
    },
]

// Mock Chart Data
const userGrowthData = [
    { name: 'Mon', users: 10 },
    { name: 'Tue', users: 15 },
    { name: 'Wed', users: 12 },
    { name: 'Thu', users: 20 },
    { name: 'Fri', users: 25 },
    { name: 'Sat', users: 30 },
    { name: 'Sun', users: 35 },
]

const topTargetsData = [
    { name: 'HER2', value: 400 },
    { name: 'LIV-1', value: 300 },
    { name: 'TROP2', value: 200 },
    { name: 'CD30', value: 100 },
]

const COLORS = ['#8b5cf6', '#ec4899', '#3b82f6', '#10b981']

// Mock Live Feed
const liveFeed = [
    { id: 1, message: 'Simulation #Job-992 failed (Timeout)', time: '2 mins ago', type: 'error' },
    { id: 2, message: 'User kim@bionet.com purchased Pro Plan', time: '15 mins ago', type: 'success' },
    { id: 3, message: 'Crawler: 150 new PubMed articles indexed', time: '1 hour ago', type: 'info' },
    { id: 4, message: 'System Backup completed successfully', time: '3 hours ago', type: 'info' },
    { id: 5, message: 'New user registered: lee@syngene.com', time: '5 hours ago', type: 'success' },
]

// System services to monitor
const systemServices = [
    { id: 'api', name: 'API Server', endpoint: `${API_BASE_URL}/health` },
    { id: 'worker', name: 'Background Worker', endpoint: `${API_BASE_URL}/api/worker/health` },
    { id: 'db', name: 'Supabase DB', endpoint: `${API_BASE_URL}/api/db/health` },
]

type ServiceStatus = 'operational' | 'degraded' | 'down' | 'checking'

function useHealthCheck() {
    const [statuses, setStatuses] = useState<Record<string, ServiceStatus>>(
        Object.fromEntries(systemServices.map(s => [s.id, 'checking']))
    )

    const checkHealth = useCallback(async () => {
        for (const service of systemServices) {
            try {
                // Mocking health check for now as endpoints might not be fully ready
                // In production, uncomment fetch
                // const response = await fetch(service.endpoint)
                // const isHealthy = response.ok
                const isHealthy = true // Mock

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
        checkHealth()
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
                <p className="text-slate-400 mt-1">System Health & Real-time Metrics</p>
            </motion.div>

            {/* System Status - Traffic Lights */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.1 }}
                className="grid grid-cols-1 md:grid-cols-3 gap-4"
            >
                {systemServices.map((service) => {
                    const status = serviceStatuses[service.id]
                    return (
                        <Card key={service.id} className="bg-slate-900 border-slate-800">
                            <CardContent className="p-4 flex items-center justify-between">
                                <div className="flex items-center gap-3">
                                    {status === 'checking' ? (
                                        <Loader2 className="w-5 h-5 text-slate-400 animate-spin" />
                                    ) : status === 'operational' ? (
                                        <div className="w-3 h-3 rounded-full bg-green-500 shadow-[0_0_10px_rgba(34,197,94,0.5)]" />
                                    ) : (
                                        <div className="w-3 h-3 rounded-full bg-red-500 shadow-[0_0_10px_rgba(239,68,68,0.5)]" />
                                    )}
                                    <span className="font-medium text-white">{service.name}</span>
                                </div>
                                <Badge variant="outline" className={`capitalize ${status === 'operational' ? 'border-green-500/30 text-green-400' :
                                        status === 'checking' ? 'border-slate-600 text-slate-400' :
                                            'border-red-500/30 text-red-400'
                                    }`}>
                                    {status}
                                </Badge>
                            </CardContent>
                        </Card>
                    )
                })}
            </motion.div>

            {/* KPI Cards */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.2 }}
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
                {/* User Growth Chart */}
                <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    transition={{ delay: 0.3 }}
                >
                    <Card className="bg-slate-900 border-slate-800 h-full">
                        <CardHeader>
                            <CardTitle className="text-white flex items-center gap-2">
                                <Users className="w-5 h-5" />
                                User Growth
                            </CardTitle>
                            <CardDescription className="text-slate-400">
                                Weekly new user registrations
                            </CardDescription>
                        </CardHeader>
                        <CardContent>
                            <div className="h-[300px] w-full">
                                <ResponsiveContainer width="100%" height="100%">
                                    <LineChart data={userGrowthData}>
                                        <CartesianGrid strokeDasharray="3 3" stroke="#334155" />
                                        <XAxis dataKey="name" stroke="#94a3b8" />
                                        <YAxis stroke="#94a3b8" />
                                        <Tooltip
                                            contentStyle={{ backgroundColor: '#1e293b', borderColor: '#334155', color: '#fff' }}
                                            itemStyle={{ color: '#fff' }}
                                        />
                                        <Line type="monotone" dataKey="users" stroke="#8b5cf6" strokeWidth={3} dot={{ r: 4, fill: '#8b5cf6' }} />
                                    </LineChart>
                                </ResponsiveContainer>
                            </div>
                        </CardContent>
                    </Card>
                </motion.div>

                {/* Top Targets Pie Chart */}
                <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    animate={{ opacity: 1, y: 0 }}
                    transition={{ delay: 0.4 }}
                >
                    <Card className="bg-slate-900 border-slate-800 h-full">
                        <CardHeader>
                            <CardTitle className="text-white flex items-center gap-2">
                                <Activity className="w-5 h-5" />
                                Top Targets
                            </CardTitle>
                            <CardDescription className="text-slate-400">
                                Most analyzed targets this month
                            </CardDescription>
                        </CardHeader>
                        <CardContent>
                            <div className="h-[300px] w-full flex items-center justify-center">
                                <ResponsiveContainer width="100%" height="100%">
                                    <PieChart>
                                        <Pie
                                            data={topTargetsData}
                                            cx="50%"
                                            cy="50%"
                                            innerRadius={60}
                                            outerRadius={100}
                                            fill="#8884d8"
                                            paddingAngle={5}
                                            dataKey="value"
                                        >
                                            {topTargetsData.map((entry, index) => (
                                                <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />
                                            ))}
                                        </Pie>
                                        <Tooltip
                                            contentStyle={{ backgroundColor: '#1e293b', borderColor: '#334155', color: '#fff' }}
                                            itemStyle={{ color: '#fff' }}
                                        />
                                    </PieChart>
                                </ResponsiveContainer>
                                <div className="space-y-2">
                                    {topTargetsData.map((entry, index) => (
                                        <div key={index} className="flex items-center gap-2">
                                            <div className="w-3 h-3 rounded-full" style={{ backgroundColor: COLORS[index % COLORS.length] }} />
                                            <span className="text-sm text-slate-300">{entry.name} ({entry.value})</span>
                                        </div>
                                    ))}
                                </div>
                            </div>
                        </CardContent>
                    </Card>
                </motion.div>
            </div>

            {/* Live Feed */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.5 }}
            >
                <Card className="bg-slate-900 border-slate-800">
                    <CardHeader>
                        <CardTitle className="text-white flex items-center gap-2">
                            <Activity className="w-5 h-5" />
                            Live Feed
                        </CardTitle>
                        <CardDescription className="text-slate-400">
                            Real-time system events and logs
                        </CardDescription>
                    </CardHeader>
                    <CardContent>
                        <div className="space-y-4">
                            {liveFeed.map((item) => (
                                <div key={item.id} className="flex items-start gap-4 p-3 rounded-lg bg-slate-800/30 border border-slate-800 hover:bg-slate-800/50 transition-colors">
                                    <div className={`mt-1 w-2 h-2 rounded-full flex-shrink-0 ${item.type === 'error' ? 'bg-red-500 shadow-[0_0_8px_rgba(239,68,68,0.5)]' :
                                            item.type === 'success' ? 'bg-green-500 shadow-[0_0_8px_rgba(34,197,94,0.5)]' :
                                                'bg-blue-500 shadow-[0_0_8px_rgba(59,130,246,0.5)]'
                                        }`} />
                                    <div className="flex-1">
                                        <p className="text-sm text-slate-200">{item.message}</p>
                                        <p className="text-xs text-slate-500 mt-1">{item.time}</p>
                                    </div>
                                </div>
                            ))}
                        </div>
                    </CardContent>
                </Card>
            </motion.div>
        </div>
    )
}
