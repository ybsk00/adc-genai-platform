import { useState, useEffect, useCallback } from 'react'
import { motion } from 'framer-motion'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import { Button } from '@/components/ui/button'
import { Progress } from '@/components/ui/progress'
import {
    Activity, Loader2, Database, Zap, Box, AlertTriangle,
    FileText, CheckCircle, Clock, ArrowRight
} from 'lucide-react'
import { API_BASE_URL } from '@/lib/api'
import { useQuery } from '@tanstack/react-query'
import { useNavigate } from 'react-router-dom'
import { getSession } from '@/lib/supabase'

// --- Types ---
interface DashboardStats {
    design_stats: {
        pending: number
        processing: number
        completed: number
    }
    data_stats: {
        golden_set_approved: number
        golden_set_candidates: number
        kb_total: number
        kb_completed: number
    }
    inventory_stats: {
        antibodies: number
        reagents: number
        targets: number
    }
    system_stats: {
        ai_cost_month: number
        data_velocity_24h: number
    }
    recent_activity: {
        type: 'paper' | 'goldenset'
        id: string
        title: string
        time: string
        desc: string
    }[]
}

export function AdminOverview() {
    const navigate = useNavigate()

    // --- Queries ---
    const { data: stats, isLoading } = useQuery({
        queryKey: ['adminStats'],
        queryFn: async () => {
            const { session } = await getSession()
            const response = await fetch(`${API_BASE_URL}/api/admin/stats`, {
                headers: { Authorization: `Bearer ${session?.access_token}` }
            })
            if (!response.ok) throw new Error('Failed to fetch stats')
            return response.json() as Promise<DashboardStats>
        },
        refetchInterval: 30000 // Refresh every 30s
    })

    // --- Render Helpers ---
    const containerVariants = {
        hidden: { opacity: 0 },
        visible: {
            opacity: 1,
            transition: {
                staggerChildren: 0.1
            }
        }
    }

    const itemVariants = {
        hidden: { y: 20, opacity: 0 },
        visible: { y: 0, opacity: 1 }
    }

    if (isLoading) {
        return <div className="flex justify-center p-20"><Loader2 className="w-8 h-8 animate-spin text-blue-500" /></div>
    }

    return (
        <motion.div
            className="space-y-6"
            variants={containerVariants}
            initial="hidden"
            animate="visible"
        >
            {/* Page Header */}
            <motion.div variants={itemVariants}>
                <h1 className="text-2xl font-bold text-white">Admin Overview</h1>
                <p className="text-slate-400 mt-1">System Health & Real-time Metrics</p>
            </motion.div>

            {/* Main Stats Grid */}
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
                {/* 1. Design Engine */}
                <motion.div variants={itemVariants}>
                    <Card className="bg-slate-900 border-slate-800 h-full">
                        <CardHeader className="pb-2">
                            <CardTitle className="text-sm font-medium text-slate-400 flex items-center gap-2">
                                <Zap className="w-4 h-4 text-yellow-500" /> Design Engine
                            </CardTitle>
                        </CardHeader>
                        <CardContent>
                            <div className="space-y-3">
                                <div className="flex justify-between items-center">
                                    <span className="text-2xl font-bold text-white">{stats?.design_stats.processing || 0}</span>
                                    <Badge variant="outline" className="text-yellow-400 border-yellow-500/30 bg-yellow-500/10 animate-pulse">
                                        Running
                                    </Badge>
                                </div>
                                <div className="grid grid-cols-2 gap-2 text-xs text-slate-500">
                                    <div className="bg-slate-950 p-2 rounded border border-slate-800">
                                        <div className="font-semibold text-slate-300">{stats?.design_stats.pending || 0}</div>
                                        <div>Pending</div>
                                    </div>
                                    <div className="bg-slate-950 p-2 rounded border border-slate-800">
                                        <div className="font-semibold text-green-400">{stats?.design_stats.completed || 0}</div>
                                        <div>Done</div>
                                    </div>
                                </div>
                            </div>
                        </CardContent>
                    </Card>
                </motion.div>

                {/* 2. Data Maturity */}
                <motion.div variants={itemVariants}>
                    <Card className="bg-slate-900 border-slate-800 h-full">
                        <CardHeader className="pb-2">
                            <CardTitle className="text-sm font-medium text-slate-400 flex items-center gap-2">
                                <Database className="w-4 h-4 text-blue-500" /> Data Maturity
                            </CardTitle>
                        </CardHeader>
                        <CardContent>
                            <div className="space-y-4">
                                <div className="flex justify-between items-center">
                                    <div>
                                        <div className="text-2xl font-bold text-white">{stats?.data_stats.golden_set_approved || 0}</div>
                                        <div className="text-xs text-slate-500">Golden Sets</div>
                                    </div>
                                    <div className="text-right cursor-pointer hover:opacity-80 transition-opacity" onClick={() => navigate('/admin/staging')}>
                                        <div className="text-xl font-bold text-blue-400">{stats?.data_stats.golden_set_candidates || 0}</div>
                                        <div className="text-xs text-blue-500 flex items-center justify-end gap-1">
                                            Candidates <ArrowRight className="w-3 h-3" />
                                        </div>
                                    </div>
                                </div>
                                <div className="h-1.5 w-full bg-slate-800 rounded-full overflow-hidden">
                                    <div
                                        className="h-full bg-blue-500"
                                        style={{ width: `${Math.min(((stats?.data_stats.kb_completed || 0) / (stats?.data_stats.kb_total || 1)) * 100, 100)}%` }}
                                    />
                                </div>
                                <div className="text-xs text-slate-500 flex justify-between">
                                    <span>KB Indexed</span>
                                    <span>{stats?.data_stats.kb_completed}/{stats?.data_stats.kb_total}</span>
                                </div>
                            </div>
                        </CardContent>
                    </Card>
                </motion.div>

                {/* 3. Inventory (with Goal Progress) */}
                <motion.div variants={itemVariants}>
                    <Card className="bg-slate-900 border-slate-800 h-full">
                        <CardHeader className="pb-2">
                            <CardTitle className="text-sm font-medium text-slate-400 flex items-center gap-2">
                                <Box className="w-4 h-4 text-purple-500" /> Inventory
                            </CardTitle>
                        </CardHeader>
                        <CardContent>
                            <div className="space-y-4">
                                <div>
                                    <div className="flex justify-between text-xs mb-1">
                                        <span className="text-slate-300">Antibodies</span>
                                        <span className="text-purple-400">{((stats?.inventory_stats.antibodies || 0) / 10000 * 100).toFixed(1)}%</span>
                                    </div>
                                    <div className="h-2 w-full bg-slate-800 rounded-full overflow-hidden">
                                        <div
                                            className="h-full bg-purple-500 transition-all duration-500"
                                            style={{ width: `${Math.min(((stats?.inventory_stats.antibodies || 0) / 10000 * 100), 100)}%` }}
                                        />
                                    </div>
                                    <div className="text-xs text-slate-500 mt-1 text-right">
                                        {stats?.inventory_stats.antibodies.toLocaleString()} / 10,000 Goal
                                    </div>
                                </div>
                                <div className="grid grid-cols-2 gap-2 text-xs">
                                    <div className="flex justify-between p-2 bg-slate-950 rounded border border-slate-800">
                                        <span className="text-slate-500">Reagents</span>
                                        <span className="text-white font-mono">{stats?.inventory_stats.reagents}</span>
                                    </div>
                                    <div className="flex justify-between p-2 bg-slate-950 rounded border border-slate-800">
                                        <span className="text-slate-500">Targets</span>
                                        <span className="text-white font-mono">{stats?.inventory_stats.targets}</span>
                                    </div>
                                </div>
                            </div>
                        </CardContent>
                    </Card>
                </motion.div>

                {/* 4. System Health & Cost */}
                <motion.div variants={itemVariants}>
                    <Card className="bg-slate-900 border-slate-800 h-full">
                        <CardHeader className="pb-2">
                            <CardTitle className="text-sm font-medium text-slate-400 flex items-center gap-2">
                                <Activity className="w-4 h-4 text-green-500" /> System Health
                            </CardTitle>
                        </CardHeader>
                        <CardContent>
                            <div className="space-y-4">
                                <div className="flex items-center justify-between p-3 bg-slate-950 rounded-lg border border-slate-800">
                                    <div className="flex items-center gap-2">
                                        <div className={`w-2 h-2 rounded-full ${(stats?.system_stats.ai_cost_month || 0) > 100 ? 'bg-red-500 animate-pulse' : 'bg-green-500'}`} />
                                        <span className="text-xs text-slate-400">AI Cost (Mo)</span>
                                    </div>
                                    <span className={`font-mono font-bold ${(stats?.system_stats.ai_cost_month || 0) > 100 ? 'text-red-500 animate-pulse' : 'text-white'}`}>
                                        ${stats?.system_stats.ai_cost_month}
                                    </span>
                                </div>
                                <div className="flex items-center justify-between p-3 bg-slate-950 rounded-lg border border-slate-800">
                                    <div className="flex items-center gap-2">
                                        <Activity className="w-3 h-3 text-blue-500" />
                                        <span className="text-xs text-slate-400">Data Velocity (24h)</span>
                                    </div>
                                    <span className="font-mono font-bold text-white">
                                        +{stats?.system_stats.data_velocity_24h} items
                                    </span>
                                </div>
                            </div>
                        </CardContent>
                    </Card>
                </motion.div>
            </div>

            {/* Recent Activity Feed */}
            <motion.div variants={itemVariants}>
                <Card className="bg-slate-900 border-slate-800">
                    <CardHeader className="pb-3 border-b border-slate-800">
                        <CardTitle className="text-sm font-medium text-slate-200 flex items-center gap-2">
                            <Clock className="w-4 h-4" /> Recent Activity (Live Feed)
                        </CardTitle>
                    </CardHeader>
                    <CardContent className="p-0">
                        <div className="divide-y divide-slate-800">
                            {stats?.recent_activity.map((item, i) => (
                                <div
                                    key={i}
                                    className="p-4 hover:bg-slate-800/50 transition-colors cursor-pointer flex items-center gap-4 group"
                                    onClick={() => {
                                        if (item.type === 'paper') navigate('/admin/knowledge')
                                        else if (item.type === 'goldenset') navigate('/admin/goldenset')
                                    }}
                                >
                                    <div className={`w-8 h-8 rounded-full flex items-center justify-center shrink-0 ${item.type === 'paper' ? 'bg-blue-900/30 text-blue-400' : 'bg-yellow-900/30 text-yellow-400'
                                        }`}>
                                        {item.type === 'paper' ? <FileText className="w-4 h-4" /> : <CheckCircle className="w-4 h-4" />}
                                    </div>
                                    <div className="flex-1 min-w-0">
                                        <p className="text-sm text-slate-200 font-medium truncate group-hover:text-blue-400 transition-colors">
                                            {item.desc}
                                        </p>
                                        <p className="text-xs text-slate-500 truncate">
                                            {item.title}
                                        </p>
                                    </div>
                                    <div className="text-xs text-slate-600 font-mono whitespace-nowrap">
                                        {new Date(item.time).toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' })}
                                    </div>
                                </div>
                            ))}
                            {(!stats?.recent_activity || stats.recent_activity.length === 0) && (
                                <div className="p-8 text-center text-slate-500 text-sm">
                                    No recent activity
                                </div>
                            )}
                        </div>
                    </CardContent>
                </Card>
            </motion.div>
        </motion.div>
    )
}

