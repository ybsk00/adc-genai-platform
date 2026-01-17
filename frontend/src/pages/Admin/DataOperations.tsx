import { useState } from 'react'
import { motion } from 'framer-motion'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import {
    Dialog,
    DialogContent,
    DialogDescription,
    DialogFooter,
    DialogHeader,
    DialogTitle,
} from '@/components/ui/dialog'
import {
    Database,
    RefreshCw,
    FileText,
    Loader2,
    CheckCircle,
    Clock,
    AlertTriangle,
    Table as TableIcon
} from 'lucide-react'
import { toast } from 'sonner'
import { GoldenSetEditor } from '@/components/admin/GoldenSetEditor'

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

type ConfirmationType = 'syncAll' | 'rebuildRAG' | null

export function DataOperations() {
    const [syncingIds, setSyncingIds] = useState<string[]>([])
    const [confirmationType, setConfirmationType] = useState<ConfirmationType>(null)
    const [isExecuting, setIsExecuting] = useState(false)
    const [progress, setProgress] = useState(0)

    /**
     * [Dev Note: Confirmation Modal]
     * Sync All / Rebuild Index는 비용+시간이 많이 듦
     * 실수로 누르지 않도록 경고 팝업 필수
     */
    const handleOpenConfirmation = (type: ConfirmationType) => {
        setConfirmationType(type)
        setProgress(0)
    }

    const handleConfirmAction = async () => {
        if (!confirmationType) return

        setIsExecuting(true)
        setProgress(0)

        try {
            if (confirmationType === 'syncAll') {
                setSyncingIds(dataSources.map(d => d.id))

                // Progress 시뮬레이션
                for (let i = 0; i <= 100; i += 10) {
                    await new Promise(r => setTimeout(r, 500))
                    setProgress(i)
                }

                toast.success('모든 데이터 소스 동기화가 완료되었습니다.')
                setSyncingIds([])

            } else if (confirmationType === 'rebuildRAG') {
                // Progress 시뮬레이션
                for (let i = 0; i <= 100; i += 5) {
                    await new Promise(r => setTimeout(r, 400))
                    setProgress(i)
                }

                toast.success('RAG 인덱스 재구축이 완료되었습니다.')
            }
        } catch (error) {
            toast.error('작업 실행에 실패했습니다.')
        } finally {
            setIsExecuting(false)
            setConfirmationType(null)
        }
    }

    const handleSync = async (sourceId: string) => {
        setSyncingIds(prev => [...prev, sourceId])

        // TODO: API 호출
        await new Promise(resolve => setTimeout(resolve, 3000))

        toast.success(`${dataSources.find(d => d.id === sourceId)?.name} 동기화 완료`)
        setSyncingIds(prev => prev.filter(id => id !== sourceId))
    }

    const getConfirmationContent = () => {
        if (confirmationType === 'syncAll') {
            return {
                title: '전체 데이터 동기화',
                description: '모든 외부 데이터 소스를 동기화합니다.',
                warning: '이 작업은 OpenAI API 비용이 발생하며, 완료까지 약 20분이 소요됩니다.',
                estimatedTime: '20분',
                estimatedCost: '~$5'
            }
        }
        return {
            title: 'RAG 인덱스 재구축',
            description: '전체 벡터 데이터베이스를 재구축합니다.',
            warning: '이 작업은 많은 OpenAI 비용이 발생하며, 완료까지 약 30분이 소요됩니다.',
            estimatedTime: '30분',
            estimatedCost: '~$15'
        }
    }

    return (
        <div className="space-y-6">
            {/* Page Header */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                className="flex items-center justify-between"
            >
                <div>
                    <h1 className="text-2xl font-bold text-white">Data Operations</h1>
                    <p className="text-slate-400 mt-1">데이터 동기화 및 RAG 파이프라인 관리</p>
                </div>
                <div className="flex gap-3">
                    <Button
                        variant="outline"
                        className="border-slate-700 text-slate-300 hover:text-white"
                        onClick={() => handleOpenConfirmation('syncAll')}
                        disabled={syncingIds.length > 0}
                    >
                        {syncingIds.length > 0 ? (
                            <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                        ) : (
                            <RefreshCw className="w-4 h-4 mr-2" />
                        )}
                        Sync All
                    </Button>
                    <Button
                        className="bg-purple-500 hover:bg-purple-600"
                        onClick={() => handleOpenConfirmation('rebuildRAG')}
                    >
                        <Database className="w-4 h-4 mr-2" />
                        Rebuild RAG Index
                    </Button>
                </div>
            </motion.div>

            <Tabs defaultValue="sources" className="space-y-4">
                <TabsList className="bg-slate-900 border border-slate-800">
                    <TabsTrigger value="sources" className="data-[state=active]:bg-slate-800 text-slate-400 data-[state=active]:text-white">
                        <Database className="w-4 h-4 mr-2" />
                        Data Sources
                    </TabsTrigger>
                    <TabsTrigger value="editor" className="data-[state=active]:bg-slate-800 text-slate-400 data-[state=active]:text-white">
                        <TableIcon className="w-4 h-4 mr-2" />
                        Golden Set Editor
                    </TabsTrigger>
                </TabsList>

                <TabsContent value="sources" className="space-y-4">
                    {/* Data Sources */}
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
                                        </CardContent>
                                    </Card>
                                </motion.div>
                            )
                        })}
                    </motion.div>

                    {/* RAG Pipeline Status */}
                    <motion.div
                        initial={{ opacity: 0, y: 20 }}
                        animate={{ opacity: 1, y: 0 }}
                        transition={{ delay: 0.3 }}
                    >
                        <Card className="bg-slate-900 border-slate-800">
                            <CardHeader>
                                <CardTitle className="text-white flex items-center gap-2">
                                    <FileText className="w-5 h-5" />
                                    RAG Pipeline Status
                                </CardTitle>
                                <CardDescription className="text-slate-400">
                                    벡터 데이터베이스 현황
                                </CardDescription>
                            </CardHeader>
                            <CardContent>
                                <div className="grid grid-cols-4 gap-6">
                                    <div className="text-center p-4 bg-slate-800/50 rounded-lg">
                                        <p className="text-3xl font-bold text-white">18,174</p>
                                        <p className="text-sm text-slate-400 mt-1">Total Chunks</p>
                                    </div>
                                    <div className="text-center p-4 bg-slate-800/50 rounded-lg">
                                        <p className="text-3xl font-bold text-white">1,536</p>
                                        <p className="text-sm text-slate-400 mt-1">Embedding Dim</p>
                                    </div>
                                    <div className="text-center p-4 bg-slate-800/50 rounded-lg">
                                        <p className="text-3xl font-bold text-white">0.92</p>
                                        <p className="text-sm text-slate-400 mt-1">Avg Relevance</p>
                                    </div>
                                    <div className="text-center p-4 bg-slate-800/50 rounded-lg">
                                        <p className="text-3xl font-bold text-green-400">Healthy</p>
                                        <p className="text-sm text-slate-400 mt-1">Index Status</p>
                                    </div>
                                </div>
                            </CardContent>
                        </Card>
                    </motion.div>
                </TabsContent>

                <TabsContent value="editor">
                    <Card className="bg-slate-900 border-slate-800">
                        <CardHeader>
                            <CardTitle className="text-white">Golden Set Editor</CardTitle>
                            <CardDescription className="text-slate-400">
                                FDA 승인 ADC 데이터를 직접 수정하고 관리합니다. (Human-in-the-Loop)
                            </CardDescription>
                        </CardHeader>
                        <CardContent>
                            <GoldenSetEditor />
                        </CardContent>
                    </Card>
                </TabsContent>
            </Tabs>

            {/* Confirmation Modal */}
            <Dialog open={!!confirmationType} onOpenChange={() => !isExecuting && setConfirmationType(null)}>
                <DialogContent className="bg-slate-900 border-slate-800 text-white">
                    <DialogHeader>
                        <DialogTitle className="flex items-center gap-2">
                            <AlertTriangle className="w-5 h-5 text-amber-500" />
                            {getConfirmationContent().title}
                        </DialogTitle>
                        <DialogDescription className="text-slate-400">
                            {getConfirmationContent().description}
                        </DialogDescription>
                    </DialogHeader>

                    {isExecuting ? (
                        <div className="py-8 space-y-4">
                            <div className="flex items-center justify-center gap-3">
                                <Loader2 className="w-6 h-6 text-purple-500 animate-spin" />
                                <span className="text-white">작업 진행 중...</span>
                            </div>
                            <div className="space-y-2">
                                <div className="flex justify-between text-sm">
                                    <span className="text-slate-400">진행률</span>
                                    <span className="text-white font-medium">{progress}%</span>
                                </div>
                                <div className="h-2 bg-slate-800 rounded-full overflow-hidden">
                                    <motion.div
                                        className="h-full bg-gradient-to-r from-purple-500 to-pink-500"
                                        initial={{ width: 0 }}
                                        animate={{ width: `${progress}%` }}
                                    />
                                </div>
                            </div>
                            <p className="text-xs text-slate-500 text-center">
                                브라우저를 닫지 마세요. 작업이 완료되면 자동으로 알림이 표시됩니다.
                            </p>
                        </div>
                    ) : (
                        <>
                            <div className="py-4 space-y-4">
                                <div className="p-4 bg-amber-500/10 border border-amber-500/30 rounded-lg">
                                    <p className="text-sm text-amber-400">
                                        ⚠️ {getConfirmationContent().warning}
                                    </p>
                                </div>
                                <div className="grid grid-cols-2 gap-4 text-sm">
                                    <div className="p-3 bg-slate-800/50 rounded-lg">
                                        <p className="text-slate-400">예상 소요시간</p>
                                        <p className="text-white font-medium text-lg">{getConfirmationContent().estimatedTime}</p>
                                    </div>
                                    <div className="p-3 bg-slate-800/50 rounded-lg">
                                        <p className="text-slate-400">예상 비용</p>
                                        <p className="text-white font-medium text-lg">{getConfirmationContent().estimatedCost}</p>
                                    </div>
                                </div>
                            </div>
                            <DialogFooter>
                                <Button
                                    variant="outline"
                                    onClick={() => setConfirmationType(null)}
                                    className="border-slate-700 text-slate-300"
                                >
                                    취소
                                </Button>
                                <Button
                                    onClick={handleConfirmAction}
                                    className="bg-amber-500 hover:bg-amber-600 text-black"
                                >
                                    확인, 실행합니다
                                </Button>
                            </DialogFooter>
                        </>
                    )}
                </DialogContent>
            </Dialog>
        </div>
    )
}
