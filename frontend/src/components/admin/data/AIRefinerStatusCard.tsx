import { useState, useEffect } from 'react'
import { motion } from 'framer-motion'
import { Card, CardContent } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Switch } from '@/components/ui/switch'
import { Progress } from '@/components/ui/progress'
import {
    Brain,
    Pause,
    Play,
    Eye,
    Loader2,
    Zap
} from 'lucide-react'
import { toast } from 'sonner'
import { API_BASE_URL } from '@/lib/api'
import {
    Dialog,
    DialogContent,
    DialogHeader,
    DialogTitle,
} from '@/components/ui/dialog'

interface DashboardData {
    system_status: 'ACTIVE' | 'PAUSED'
    cost: {
        daily_usage_usd: number
        daily_limit_usd: number
        remaining_usd: number
        is_over_limit: boolean
    }
    queue: {
        pending: number
        enriched_today: number
    }
}

interface RecentLog {
    id: string
    name: string
    outcome_type: string
    failure_reason: string | null
    properties: Record<string, unknown>
    smiles_code: string | null
    created_at: string
}

import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs"

// ... (imports remain same)

export function AIRefinerStatusCard() {
    const [loading, setLoading] = useState(true)
    const [toggling, setToggling] = useState(false)
    const [running, setRunning] = useState(false) // Clinical Refiner 상태
    const [knowledgeRunning, setKnowledgeRunning] = useState(false) // Knowledge Refiner 상태
    const [dashboard, setDashboard] = useState<DashboardData | null>(null)
    const [spotCheckOpen, setSpotCheckOpen] = useState(false)
    const [recentLogs, setRecentLogs] = useState<RecentLog[]>([])
    const [logsLoading, setLogsLoading] = useState(false)

    const fetchDashboard = async () => {
        try {
            const res = await fetch(`${API_BASE_URL}/api/admin/refiner/dashboard`)
            if (res.ok) {
                const data = await res.json()
                setDashboard(data)
            }
        } catch (error) {
            console.error('Failed to fetch dashboard', error)
        } finally {
            setLoading(false)
        }
    }

    const toggleSystem = async () => {
        setToggling(true)
        try {
            const action = dashboard?.system_status === 'ACTIVE' ? 'pause' : 'resume'
            const res = await fetch(`${API_BASE_URL}/api/admin/refiner/system/${action}`, { method: 'POST' })
            if (res.ok) {
                toast.success(`System ${action}d successfully`)
                fetchDashboard()
            }
        } catch (error) {
            toast.error('Failed to toggle system')
        } finally {
            setToggling(false)
        }
    }

    const runRefinerNow = async () => {
        setRunning(true)
        try {
            const res = await fetch(`${API_BASE_URL}/api/admin/refiner/run?limit=50`, { method: 'POST' })
            if (res.ok) {
                toast.success('AI Refiner triggered successfully')
                fetchDashboard()
            }
        } catch (error) {
            toast.error('Failed to trigger refiner')
        } finally {
            setRunning(false)
        }
    }

    const fetchRecentLogs = async (type: 'clinical' | 'pubmed' = 'clinical') => {
        setLogsLoading(true)
        try {
            // endpoint needs to support type filtering
            const endpoint = type === 'clinical'
                ? `${API_BASE_URL}/api/admin/refiner/logs`
                : `${API_BASE_URL}/api/admin/knowledge/logs` // Assuming this endpoint exists or will be created

            const res = await fetch(endpoint)
            if (res.ok) {
                const data = await res.json()
                setRecentLogs(data)
            }
        } catch (error) {
            console.error('Failed to fetch logs', error)
        } finally {
            setLogsLoading(false)
        }
    }

    const openSpotCheck = (type: 'clinical' | 'pubmed') => {
        setSpotCheckOpen(true)
        fetchRecentLogs(type)
    }

    useEffect(() => {
        fetchDashboard()
        const interval = setInterval(fetchDashboard, 30000)
        return () => clearInterval(interval)
    }, [])

    const runKnowledgeRefiner = async () => {
        setKnowledgeRunning(true)
        try {
            const res = await fetch(`${API_BASE_URL}/api/admin/knowledge/refine?limit=50`, {
                method: 'POST'
            })
            if (res.ok) {
                toast.success('📚 Knowledge Refiner가 실행되었습니다! (50건)')
            } else {
                toast.error('실행 실패')
            }
        } catch {
            toast.error('네트워크 오류')
        } finally {
            setKnowledgeRunning(false)
        }
    }

    if (loading) {
        return (
            <Card className="bg-gradient-to-br from-purple-900/50 to-slate-900 border-purple-500/30">
                <CardContent className="p-6 flex items-center justify-center min-h-[200px]">
                    <Loader2 className="w-6 h-6 animate-spin text-purple-400" />
                </CardContent>
            </Card>
        )
    }

    if (!dashboard) return null

    const costPercent = (dashboard.cost.daily_usage_usd / dashboard.cost.daily_limit_usd) * 100
    const costColor = costPercent >= 80 ? 'bg-red-500' : costPercent >= 50 ? 'bg-yellow-500' : 'bg-green-500'

    return (
        <>
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.05 }}
            >
                <Card className="bg-gradient-to-br from-purple-900/50 to-slate-900 border-purple-500/30">
                    <CardContent className="p-6">
                        {/* Header */}
                        <div className="flex items-center justify-between mb-4">
                            <h3 className="text-lg font-semibold text-white flex items-center gap-2">
                                <Brain className="w-5 h-5 text-purple-400" />
                                AI Refiner Status
                                {dashboard.system_status === 'ACTIVE' && (
                                    <span className="w-2 h-2 bg-green-500 rounded-full animate-pulse" />
                                )}
                            </h3>
                            <div className="flex items-center gap-2">
                                {dashboard.system_status === 'PAUSED' ? (
                                    <Pause className="w-4 h-4 text-yellow-400" />
                                ) : (
                                    <Play className="w-4 h-4 text-green-400" />
                                )}
                                <Switch
                                    checked={dashboard.system_status === 'ACTIVE'}
                                    onCheckedChange={toggleSystem}
                                    disabled={toggling}
                                />
                            </div>
                        </div>

                        <Tabs defaultValue="clinical" className="w-full">
                            <TabsList className="grid w-full grid-cols-2 mb-4 bg-slate-800/50">
                                <TabsTrigger value="clinical">Clinical Trials</TabsTrigger>
                                <TabsTrigger value="pubmed">PubMed / BioRxiv</TabsTrigger>
                            </TabsList>

                            <TabsContent value="clinical" className="space-y-4">
                                {/* Cost Progress */}
                                <div>
                                    <div className="flex justify-between text-sm mb-1">
                                        <span className="text-slate-400">Today's Cost</span>
                                        <span className={costPercent >= 80 ? 'text-red-400' : 'text-slate-300'}>
                                            ${dashboard.cost.daily_usage_usd.toFixed(2)} / ${dashboard.cost.daily_limit_usd}
                                        </span>
                                    </div>
                                    <Progress
                                        value={costPercent}
                                        className="h-2"
                                    />
                                    <div className={`h-2 rounded-full ${costColor}`} style={{ width: `${Math.min(costPercent, 100)}%`, marginTop: '-8px' }} />
                                </div>

                                {/* Queue Stats */}
                                <div className="grid grid-cols-2 gap-4 text-center">
                                    <div className="bg-slate-800/50 rounded-lg p-3">
                                        <p className="text-2xl font-bold text-orange-400">{dashboard.queue.pending.toLocaleString()}</p>
                                        <p className="text-xs text-slate-400">Pending</p>
                                    </div>
                                    <div className="bg-slate-800/50 rounded-lg p-3">
                                        <p className="text-2xl font-bold text-green-400">{dashboard.queue.enriched_today.toLocaleString()}</p>
                                        <p className="text-xs text-slate-400">Enriched Today</p>
                                    </div>
                                </div>

                                {/* Action Buttons */}
                                <div className="grid grid-cols-2 gap-3">
                                    <Button
                                        variant="outline"
                                        className="border-purple-500/30 text-purple-300 hover:bg-purple-500/10"
                                        onClick={runRefinerNow}
                                        disabled={running || dashboard.system_status === 'PAUSED'}
                                    >
                                        {running ? (
                                            <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                                        ) : (
                                            <Zap className="w-4 h-4 mr-2" />
                                        )}
                                        Run Now
                                    </Button>
                                    <Button
                                        variant="outline"
                                        className="w-full border-purple-500/30 text-purple-300 hover:bg-purple-500/10"
                                        onClick={openSpotCheck}
                                    >
                                        <Eye className="w-4 h-4 mr-2" />
                                        Spot Check
                                    </Button>
                                </div>
                            </TabsContent>

                            <TabsContent value="pubmed" className="space-y-4">
                                <div className="bg-slate-800/50 rounded-lg p-4 text-center border border-blue-500/20">
                                    <p className="text-sm text-blue-300 font-medium mb-1">Knowledge Base Refiner</p>
                                    <p className="text-xs text-slate-400">
                                        Analyzes PubMed/BioRxiv abstracts using Gemini Flash.<br />
                                        Extracts Summary, Relevance Score, and AI Reasoning.
                                    </p>
                                </div>

                                <div className="grid grid-cols-2 gap-3">
                                    <Button
                                        variant="outline"
                                        className="border-blue-500/30 text-blue-300 hover:bg-blue-500/10"
                                        onClick={runKnowledgeRefiner}
                                        disabled={knowledgeRunning}
                                    >
                                        {knowledgeRunning ? (
                                            <Loader2 className="w-5 h-5 mr-2 animate-spin" />
                                        ) : (
                                            <Zap className="w-5 h-5 mr-2" />
                                        )}
                                        Run Now
                                    </Button>
                                    <Button
                                        variant="outline"
                                        className="border-blue-500/30 text-blue-300 hover:bg-blue-500/10"
                                        onClick={() => openSpotCheck('pubmed')}
                                    >
                                        <Eye className="w-4 h-4 mr-2" />
                                        Spot Check
                                    </Button>
                                </div>

                                <div className="text-center mt-2">
                                    <p className="text-xs text-slate-500">
                                        * Automatically runs daily via scheduler.<br />
                                        * Use this button for immediate processing of pending items.
                                    </p>
                                </div>
                            </TabsContent>
                        </Tabs>
                    </CardContent>
                </Card>
            </motion.div>

            {/* Spot Check Modal (Existing) */}
            <Dialog open={spotCheckOpen} onOpenChange={setSpotCheckOpen}>
                <DialogContent className="max-w-2xl bg-slate-900 border-slate-700">
                    <DialogHeader>
                        <DialogTitle className="text-white">🕵️ Recent AI Refiner Results</DialogTitle>
                    </DialogHeader>

                    {logsLoading ? (
                        <div className="flex justify-center py-8">
                            <Loader2 className="w-6 h-6 animate-spin text-purple-400" />
                        </div>
                    ) : (
                        <div className="space-y-3 max-h-96 overflow-y-auto">
                            {recentLogs.map((log) => (
                                <div key={log.id} className="bg-slate-800 rounded-lg p-4">
                                    <div className="flex justify-between items-start mb-2">
                                        <h4 className="font-medium text-white">{log.name || 'Unknown'}</h4>
                                        <span className={`text-xs px-2 py-1 rounded ${log.outcome_type === 'Success' ? 'bg-green-500/20 text-green-400' :
                                            log.outcome_type === 'Failure' ? 'bg-red-500/20 text-red-400' :
                                                'bg-slate-500/20 text-slate-400'
                                            }`}>
                                            {log.outcome_type || 'Unknown'}
                                        </span>
                                    </div>
                                    {log.failure_reason && (
                                        <p className="text-sm text-red-400 mb-2">Reason: {log.failure_reason}</p>
                                    )}
                                    {log.smiles_code && (
                                        <p className="text-xs text-slate-400 font-mono truncate">SMILES: {log.smiles_code}</p>
                                    )}
                                    <p className="text-xs text-slate-500 mt-2">
                                        {new Date(log.created_at).toLocaleString()}
                                    </p>
                                </div>
                            ))}
                            {recentLogs.length === 0 && (
                                <p className="text-center text-slate-400 py-4">No processed records yet.</p>
                            )}
                        </div>
                    )}
                </DialogContent>
            </Dialog>
        </>
    )
}
