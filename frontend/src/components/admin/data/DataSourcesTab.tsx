import { useState, useRef } from 'react'
import { motion } from 'framer-motion'
import { Card, CardContent } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import {
    RefreshCw,
    Loader2,
    CheckCircle,
    Clock,
    AlertTriangle,
    Settings,
    Square,
    ShieldAlert
} from 'lucide-react'
import { toast } from 'sonner'
import { DataSourceSettingsDialog } from './DataSourceSettingsDialog'
import { AIRefinerStatusCard } from './AIRefinerStatusCard'

interface DataSource {
    id: string
    name: string
    description: string
    lastSync: string
    status: 'synced' | 'syncing' | 'error'
    recordCount: number
    estimatedTime: string
    endpoint?: string
    jobId?: string
}

import { API_BASE_URL } from '@/lib/api'

export function DataSourcesTab() {
    const [syncingIds, setSyncingIds] = useState<string[]>([])
    const [settingsOpen, setSettingsOpen] = useState(false)
    const [selectedSource, setSelectedSource] = useState<{ id: string, name: string } | null>(null)

    // Polling intervals ref for cleanup
    const pollingIntervalsRef = useRef<Map<string, ReturnType<typeof setInterval>>>(new Map())
    const [sources, setSources] = useState<DataSource[]>([
        {
            id: 'bulk_import',
            name: 'ClinicalTrials.gov (API v2)',
            description: 'Sync latest ADC trials via API v2',
            lastSync: 'Not synced',
            status: 'synced',
            recordCount: 0,
            estimatedTime: '1-3 min',
            endpoint: `${API_BASE_URL}/api/scheduler/bulk/import`
        },
        {
            id: 'pubmed',
            name: 'PubMed/BioRxiv',
            description: 'Scientific literature → knowledge_base',
            lastSync: 'Not synced',
            status: 'synced',
            recordCount: 0,
            estimatedTime: '15 min',
            endpoint: `${API_BASE_URL}/api/scheduler/sync/pubmed`
        },
        {
            id: 'openfda',
            name: 'OpenFDA Approved Labels',
            description: 'FDA-approved ADCs → golden_set_library',
            lastSync: 'Not synced',
            status: 'synced',
            recordCount: 0,
            estimatedTime: '3 min',
            endpoint: `${API_BASE_URL}/api/scheduler/sync/openfda`
        },
        {
            id: 'creative',
            name: 'Creative Biolabs',
            description: 'Commercial reagents → commercial_reagents',
            lastSync: 'Not synced',
            status: 'synced',
            recordCount: 0,
            estimatedTime: '10 min',
            endpoint: `${API_BASE_URL}/api/scheduler/crawler/creative/run`
        },
        {
            id: 'ai_refiner',
            name: 'AI Refiner',
            description: 'Process pending records with LLM analysis',
            lastSync: 'Not synced',
            status: 'synced',
            recordCount: 0,
            estimatedTime: '5 min',
            endpoint: `${API_BASE_URL}/api/scheduler/refiner/run`
        },
    ])



    const pollJobStatus = async (jobId: string, sourceId: string) => {
        // Clear any existing interval for this source
        const existingInterval = pollingIntervalsRef.current.get(sourceId)
        if (existingInterval) {
            clearInterval(existingInterval)
        }

        const interval = setInterval(async () => {
            try {
                const res = await fetch(`${API_BASE_URL}/api/scheduler/sync/${jobId}`)
                if (res.ok) {
                    const status = await res.json()
                    setSources(prev => prev.map(s => s.id === sourceId ? {
                        ...s,
                        status: status.status === 'running' ? 'syncing' :
                            status.status === 'completed' ? 'synced' :
                                status.status === 'failed' ? 'error' :
                                    status.status === 'stopped' ? 'synced' : s.status,
                        recordCount: status.records_drafted || s.recordCount,
                        lastSync: status.completed_at ? 'Just now' : s.lastSync,
                        jobId: jobId
                    } : s))

                    if (status.status === 'completed' || status.status === 'failed' || status.status === 'stopped') {
                        clearInterval(interval)
                        pollingIntervalsRef.current.delete(sourceId)
                        setSyncingIds(prev => prev.filter(id => id !== sourceId))
                    }
                }
            } catch (e) {
                console.error("Polling failed", e)
                clearInterval(interval)
                pollingIntervalsRef.current.delete(sourceId)
                setSyncingIds(prev => prev.filter(id => id !== sourceId))
            }
        }, 2000)

        // Store interval reference for cleanup
        pollingIntervalsRef.current.set(sourceId, interval)
    }

    const fetchRecordCounts = async () => {
        try {
            const res = await fetch(`${API_BASE_URL}/api/library/stats`)
            if (res.ok) {
                const stats = await res.json()
                setSources(prev => prev.map(s => ({
                    ...s,
                    recordCount: stats[s.id] || s.recordCount,
                    lastSync: stats[`${s.id}_last_sync`] || s.lastSync
                })))
            }
        } catch (e) {
            console.error("Failed to fetch stats", e)
        }
    }

    const handleSync = async (sourceId: string, endpoint?: string, params?: Record<string, string>) => {
        if (!endpoint) {
            toast.info('This source is coming soon.')
            return
        }

        setSyncingIds(prev => [...prev, sourceId])

        try {
            const url = new URL(endpoint)
            if (params) {
                Object.entries(params).forEach(([key, value]) => {
                    url.searchParams.append(key, value)
                })
            }

            const response = await fetch(url.toString(), { method: 'POST' })
            if (!response.ok) throw new Error('Sync failed')

            const result = await response.json()
            toast.success(`${result.message}`)

            // Start polling
            if (result.job_id) {
                setSources(prev => prev.map(s => s.id === sourceId ? { ...s, jobId: result.job_id } : s))
                pollJobStatus(result.job_id, sourceId)
            } else {
                setTimeout(() => {
                    setSyncingIds(prev => prev.filter(id => id !== sourceId))
                }, 2000)
            }

        } catch (error) {
            console.error(error)
            toast.error('Failed to start sync')
            setSyncingIds(prev => prev.filter(id => id !== sourceId))
        }
    }

    const handleStop = async (sourceId: string, jobId?: string) => {
        if (!jobId) return

        try {
            const res = await fetch(`${API_BASE_URL}/api/scheduler/sync/${jobId}/stop`, { method: 'POST' })
            if (res.ok) {
                toast.info('Stopping worker...')
            }
        } catch (e) {
            toast.error('Failed to stop worker')
        }
    }

    const handleResetAll = async () => {
        if (!confirm("Are you sure you want to force stop ALL workers and clear locks? This should only be used if workers are stuck.")) return

        try {
            const res = await fetch(`${API_BASE_URL}/api/scheduler/workers/reset`, { method: 'POST' })
            if (res.ok) {
                toast.success('All workers reset and locks cleared.')

                // Clear all polling intervals
                pollingIntervalsRef.current.forEach((interval, sourceId) => {
                    clearInterval(interval)
                    console.log(`Cleared polling for ${sourceId}`)
                })
                pollingIntervalsRef.current.clear()

                // Reset all syncing states
                setSyncingIds([])

                // Reset all source statuses to synced
                setSources(prev => prev.map(s => ({
                    ...s,
                    status: 'synced' as const,
                    jobId: undefined
                })))

                // Refresh counts
                fetchRecordCounts()
            } else {
                throw new Error('Reset failed')
            }
        } catch (e) {
            toast.error('Failed to reset workers')
        }
    }

    const openSettings = (source: DataSource) => {
        setSelectedSource({ id: source.id, name: source.name })
        setSettingsOpen(true)
    }

    return (
        <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ delay: 0.1 }}
            className="grid gap-4"
        >
            {/* AI Refiner Status Card */}
            <div className="flex justify-between items-center">
                <h2 className="text-xl font-bold text-white">Data Sources</h2>
                <Button
                    variant="destructive"
                    size="sm"
                    onClick={handleResetAll}
                    className="bg-red-500/10 text-red-500 hover:bg-red-500/20 border border-red-500/20"
                >
                    <ShieldAlert className="w-4 h-4 mr-2" />
                    Force Stop All Workers
                </Button>
            </div>
            <AIRefinerStatusCard />

            {sources.map((source, index) => {
                const isSyncing = syncingIds.includes(source.id) || source.status === 'syncing'

                return (
                    <motion.div
                        key={source.id}
                        initial={{ opacity: 0, x: -20 }}
                        animate={{ opacity: 1, x: 0 }}
                        transition={{ delay: 0.1 + index * 0.05 }}
                    >
                        <Card className="bg-slate-900 border-slate-800">
                            <CardContent className="p-6">
                                <div className="flex items-center justify-between">
                                    <div className="flex items-center gap-4">
                                        <div className={`w-12 h-12 rounded-xl flex items-center justify-center ${source.status === 'error' ? 'bg-red-500/20' :
                                            isSyncing ? 'bg-blue-500/20' :
                                                'bg-green-500/20'
                                            }`}>
                                            {source.status === 'error' ? (
                                                <AlertTriangle className="w-6 h-6 text-red-400" />
                                            ) : isSyncing ? (
                                                <Loader2 className="w-6 h-6 text-blue-400 animate-spin" />
                                            ) : (
                                                <CheckCircle className="w-6 h-6 text-green-400" />
                                            )}
                                        </div>
                                        <div>
                                            <h3 className="text-lg font-semibold text-white">{source.name}</h3>
                                            <p className="text-sm text-slate-400">{source.description}</p>
                                        </div>
                                    </div>

                                    <div className="flex items-center gap-8">
                                        <div className="text-right">
                                            <p className="text-2xl font-bold text-white">{source.recordCount.toLocaleString()}</p>
                                            <p className="text-xs text-slate-400">records</p>
                                        </div>

                                        <div className="text-right min-w-[140px]">
                                            <div className="flex items-center gap-1 justify-end text-slate-400">
                                                <Clock className="w-3 h-3" />
                                                <span className="text-xs">{source.lastSync}</span>
                                            </div>
                                            <Badge
                                                variant="outline"
                                                className={`mt-1 ${source.status === 'error' ? 'border-red-500/30 text-red-400' :
                                                    isSyncing ? 'border-blue-500/30 text-blue-400' :
                                                        'border-green-500/30 text-green-400'
                                                    }`}
                                            >
                                                {isSyncing ? 'Syncing...' : source.status}
                                            </Badge>
                                        </div>

                                        <div className="flex gap-2">
                                            <Button
                                                variant="outline"
                                                size="sm"
                                                className="border-slate-700 text-slate-300 hover:text-white"
                                                onClick={() => openSettings(source)}
                                            >
                                                <Settings className="w-4 h-4" />
                                            </Button>

                                            {source.id === 'bulk_import' || source.id === 'pubmed' || source.id === 'openfda' ? (
                                                <>
                                                    <Button
                                                        variant="outline"
                                                        size="sm"
                                                        className="border-slate-700 text-slate-300 hover:text-white"
                                                        onClick={() => isSyncing ? handleStop(source.id, source.jobId) : handleSync(source.id, source.endpoint, { mode: 'daily' })}
                                                        disabled={isSyncing}
                                                    >
                                                        {isSyncing ? <Loader2 className="w-4 h-4 animate-spin" /> : "Daily"}
                                                    </Button>
                                                    <Button
                                                        variant="outline"
                                                        size="sm"
                                                        className="border-slate-700 text-slate-300 hover:text-white"
                                                        onClick={() => isSyncing ? handleStop(source.id, source.jobId) : handleSync(source.id, source.endpoint, { mode: 'full' })}
                                                        disabled={isSyncing}
                                                    >
                                                        {isSyncing ? <Loader2 className="w-4 h-4 animate-spin" /> : "Full Load"}
                                                    </Button>
                                                </>
                                            ) : (
                                                <Button
                                                    variant="outline"
                                                    size="sm"
                                                    className="border-slate-700 text-slate-300 hover:text-white"
                                                    onClick={() => isSyncing ? handleStop(source.id, source.jobId) : handleSync(source.id, source.endpoint)}
                                                >
                                                    {isSyncing ? (
                                                        <Square className="w-4 h-4 text-red-400 fill-red-400" />
                                                    ) : (
                                                        <RefreshCw className="w-4 h-4" />
                                                    )}
                                                </Button>
                                            )}
                                        </div>
                                    </div>
                                </div>
                            </CardContent>
                        </Card>
                    </motion.div>
                )
            })}

            {selectedSource && (
                <DataSourceSettingsDialog
                    open={settingsOpen}
                    onOpenChange={setSettingsOpen}
                    sourceId={selectedSource.id}
                    sourceName={selectedSource.name}
                />
            )}
        </motion.div>
    )
}
