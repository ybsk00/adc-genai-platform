import { useState } from 'react'
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
    Settings
} from 'lucide-react'
import { toast } from 'sonner'

interface DataSource {
    id: string
    name: string
    description: string
    lastSync: string
    status: 'synced' | 'syncing' | 'error'
    recordCount: number
    estimatedTime: string
}

const dataSources: DataSource[] = [
    {
        id: 'clinical',
        name: 'ClinicalTrials.gov',
        description: '임상시험 데이터',
        lastSync: '2026-01-17 03:00',
        status: 'synced',
        recordCount: 2847,
        estimatedTime: '5분'
    },
    {
        id: 'pubmed',
        name: 'PubMed/BioRxiv',
        description: '논문 데이터',
        lastSync: '2026-01-17 02:00',
        status: 'synced',
        recordCount: 15420,
        estimatedTime: '15분'
    },
    {
        id: 'news',
        name: 'Perplexity News',
        description: 'ADC 관련 뉴스',
        lastSync: '2026-01-17 08:00',
        status: 'synced',
        recordCount: 892,
        estimatedTime: '3분'
    },
    {
        id: 'goldenset',
        name: 'Golden Set Library',
        description: 'FDA 승인 ADC 데이터',
        lastSync: '2026-01-15 10:00',
        status: 'synced',
        recordCount: 15,
        estimatedTime: '1분'
    },
]

export function DataSourcesTab() {
    const [syncingIds, setSyncingIds] = useState<string[]>([])

    const handleSync = async (sourceId: string) => {
        setSyncingIds(prev => [...prev, sourceId])

        // TODO: API 호출
        await new Promise(resolve => setTimeout(resolve, 3000))

        toast.success(`${dataSources.find(d => d.id === sourceId)?.name} 동기화 완료`)
        setSyncingIds(prev => prev.filter(id => id !== sourceId))
    }

    return (
        <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ delay: 0.1 }}
            className="grid gap-4"
        >
            {dataSources.map((source, index) => {
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
                                            >
                                                <Settings className="w-4 h-4" />
                                            </Button>
                                            <Button
                                                variant="outline"
                                                size="sm"
                                                className="border-slate-700 text-slate-300 hover:text-white"
                                                onClick={() => handleSync(source.id)}
                                                disabled={isSyncing}
                                            >
                                                {isSyncing ? (
                                                    <Loader2 className="w-4 h-4 animate-spin" />
                                                ) : (
                                                    <RefreshCw className="w-4 h-4" />
                                                )}
                                            </Button>
                                        </div>
                                    </div>
                                </div>
                            </CardContent>
                        </Card>
                    </motion.div>
                )
            })}
        </motion.div>
    )
}
