import { useState, useEffect } from 'react'
import {
    Card,
    CardContent,
    CardDescription,
    CardHeader,
    CardTitle,
} from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import {
    Tooltip,
    TooltipContent,
    TooltipProvider,
    TooltipTrigger,
} from '@/components/ui/tooltip'
import {
    Table,
    TableBody,
    TableCell,
    TableHead,
    TableHeader,
    TableRow,
} from '@/components/ui/table'
import {
    FileText,
    CheckCircle,
    Trash2,
    Zap,
    ExternalLink,
    Loader2,
    Info
} from 'lucide-react'
import { toast } from 'sonner'
import { format } from 'date-fns'

interface KnowledgeBaseItem {
    id: number
    source_type: string
    title: string
    summary: string
    relevance_score: number
    source_tier: number
    rag_status: 'pending' | 'indexed' | 'excluded'
    created_at: string
    url?: string
    reasoning?: string // AI reasoning for the score
}

export function KnowledgeBaseTab() {
    const [items, setItems] = useState<KnowledgeBaseItem[]>([])
    const [loading, setLoading] = useState(true)

    const fetchData = async () => {
        setLoading(true)
        try {
            // Mock Data for UI development
            // In production: const response = await fetch('/api/knowledge-base')
            await new Promise(r => setTimeout(r, 1000))

            const mockData: KnowledgeBaseItem[] = [
                {
                    id: 1,
                    source_type: 'PubMed',
                    title: 'Efficacy of T-DXd in HER2-low metastatic breast cancer',
                    summary: 'Trastuzumab deruxtecan (T-DXd) showed significant improvement in progression-free survival compared to chemotherapy...',
                    relevance_score: 95,
                    source_tier: 1,
                    rag_status: 'pending',
                    created_at: '2026-01-18T09:00:00Z',
                    url: 'https://pubmed.ncbi.nlm.nih.gov/12345678',
                    reasoning: 'Directly discusses clinical efficacy of a major ADC drug (Enhertu) in a key indication.'
                },
                {
                    id: 2,
                    source_type: 'News',
                    title: 'Competitor X fails Phase 2 trial for TROP2 ADC',
                    summary: 'Company X announced that their TROP2-directed ADC failed to meet primary endpoints in the LUNG-01 study due to toxicity...',
                    relevance_score: 88,
                    source_tier: 2,
                    rag_status: 'pending',
                    created_at: '2026-01-18T10:30:00Z',
                    url: 'https://biopharma-news.com/article/123',
                    reasoning: 'Critical market intelligence regarding TROP2 target validation and toxicity risks.'
                },
                {
                    id: 3,
                    source_type: 'Report',
                    title: 'Global ADC Market Forecast 2026-2030',
                    summary: 'The ADC market is projected to grow at a CAGR of 15%, driven by new approvals in solid tumors...',
                    relevance_score: 60,
                    source_tier: 3,
                    rag_status: 'excluded',
                    created_at: '2026-01-17T15:00:00Z',
                    reasoning: 'General market data, less relevant for specific drug discovery tasks.'
                }
            ]
            setItems(mockData)
        } catch (error) {
            console.error(error)
            toast.error('Failed to load knowledge base')
        } finally {
            setLoading(false)
        }
    }

    useEffect(() => {
        fetchData()
    }, [])

    const handleApprove = async (id: number) => {
        // TODO: API Call
        toast.success('Item approved for indexing')
        setItems(prev => prev.map(item => item.id === id ? { ...item, rag_status: 'indexed' } : item))
    }

    const handleTrash = async (id: number) => {
        // TODO: API Call
        toast.success('Item moved to trash')
        setItems(prev => prev.filter(item => item.id !== id))
    }

    const handleAutoApprove = async () => {
        // TODO: API Call
        toast.success('Auto-approved 5 Tier-1 items')
        fetchData()
    }

    return (
        <div className="space-y-4">
            {/* Header Actions */}
            <div className="flex justify-between items-center bg-slate-900 p-4 rounded-lg border border-slate-800">
                <div>
                    <h3 className="text-lg font-medium text-white">Pending Review</h3>
                    <p className="text-sm text-slate-400">
                        {items.filter(i => i.rag_status === 'pending').length} items waiting for curation
                    </p>
                </div>
                <Button
                    className="bg-purple-500 hover:bg-purple-600 text-white"
                    onClick={handleAutoApprove}
                >
                    <Zap className="w-4 h-4 mr-2 fill-current" />
                    Auto-Approve Tier 1
                </Button>
            </div>

            {/* List View */}
            <div className="rounded-md border border-slate-800 bg-slate-900">
                <Table>
                    <TableHeader>
                        <TableRow className="border-slate-800 hover:bg-transparent">
                            <TableHead className="text-slate-400">Source</TableHead>
                            <TableHead className="text-slate-400 w-[40%]">Title & Summary</TableHead>
                            <TableHead className="text-slate-400">Relevance</TableHead>
                            <TableHead className="text-slate-400">Tier</TableHead>
                            <TableHead className="text-slate-400">Date</TableHead>
                            <TableHead className="text-slate-400 text-right">Actions</TableHead>
                        </TableRow>
                    </TableHeader>
                    <TableBody>
                        {loading ? (
                            <TableRow>
                                <TableCell colSpan={6} className="h-24 text-center">
                                    <div className="flex justify-center items-center text-slate-500">
                                        <Loader2 className="w-6 h-6 animate-spin mr-2" />
                                        Loading knowledge base...
                                    </div>
                                </TableCell>
                            </TableRow>
                        ) : items.length > 0 ? (
                            items.map((item) => (
                                <TableRow key={item.id} className="border-slate-800 hover:bg-slate-800/50">
                                    <TableCell>
                                        <Badge variant="outline" className="border-slate-700 text-slate-400">
                                            {item.source_type}
                                        </Badge>
                                    </TableCell>
                                    <TableCell>
                                        <div className="space-y-1">
                                            <div className="flex items-center gap-2">
                                                <span className="font-medium text-white line-clamp-1" title={item.title}>
                                                    {item.title}
                                                </span>
                                                {item.url && (
                                                    <a href={item.url} target="_blank" rel="noreferrer" className="text-slate-500 hover:text-purple-400">
                                                        <ExternalLink className="w-3 h-3" />
                                                    </a>
                                                )}
                                            </div>
                                            <p className="text-xs text-slate-500 line-clamp-2">
                                                {item.summary}
                                            </p>
                                        </div>
                                    </TableCell>
                                    <TableCell>
                                        <div className="flex items-center gap-2">
                                            <span className={`font-bold ${item.relevance_score >= 90 ? 'text-green-400' :
                                                    item.relevance_score >= 70 ? 'text-blue-400' :
                                                        'text-slate-400'
                                                }`}>
                                                {item.relevance_score}
                                            </span>
                                            {item.reasoning && (
                                                <TooltipProvider>
                                                    <Tooltip>
                                                        <TooltipTrigger>
                                                            <Info className="w-3 h-3 text-slate-600 hover:text-slate-400" />
                                                        </TooltipTrigger>
                                                        <TooltipContent className="bg-slate-800 border-slate-700 text-slate-200 max-w-xs">
                                                            <p className="text-xs font-semibold mb-1">AI Reasoning:</p>
                                                            <p className="text-xs">{item.reasoning}</p>
                                                        </TooltipContent>
                                                    </Tooltip>
                                                </TooltipProvider>
                                            )}
                                        </div>
                                    </TableCell>
                                    <TableCell>
                                        <div className="flex gap-1">
                                            {[1, 2, 3].map((tier) => (
                                                <div
                                                    key={tier}
                                                    className={`w-2 h-2 rounded-full ${tier <= (4 - item.source_tier) // Tier 1 = 3 dots, Tier 3 = 1 dot
                                                            ? 'bg-purple-500'
                                                            : 'bg-slate-700'
                                                        }`}
                                                />
                                            ))}
                                        </div>
                                        <span className="text-xs text-slate-500 mt-1 block">Tier {item.source_tier}</span>
                                    </TableCell>
                                    <TableCell>
                                        <span className="text-xs text-slate-500">
                                            {format(new Date(item.created_at), 'MMM d, HH:mm')}
                                        </span>
                                    </TableCell>
                                    <TableCell className="text-right">
                                        <div className="flex justify-end gap-2">
                                            <Button
                                                variant="ghost"
                                                size="icon"
                                                className="text-slate-400 hover:text-red-400 hover:bg-red-500/10"
                                                onClick={() => handleTrash(item.id)}
                                            >
                                                <Trash2 className="w-4 h-4" />
                                            </Button>
                                            <Button
                                                variant="ghost"
                                                size="icon"
                                                className="text-green-500 hover:text-green-400 hover:bg-green-500/10"
                                                onClick={() => handleApprove(item.id)}
                                            >
                                                <CheckCircle className="w-4 h-4" />
                                            </Button>
                                        </div>
                                    </TableCell>
                                </TableRow>
                            ))
                        ) : (
                            <TableRow>
                                <TableCell colSpan={6} className="h-24 text-center text-slate-500">
                                    No items found.
                                </TableCell>
                            </TableRow>
                        )}
                    </TableBody>
                </Table>
            </div>
        </div>
    )
}
