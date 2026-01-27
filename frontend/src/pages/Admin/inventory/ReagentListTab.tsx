import { useState, useEffect } from 'react'
import {
    Table,
    TableBody,
    TableCell,
    TableHead,
    TableHeader,
    TableRow
} from '@/components/ui/table'
import { Input } from '@/components/ui/input'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import { Search, Loader2, ChevronLeft, ChevronRight, AlertTriangle, CheckCircle2, Sparkles, UserCheck, Filter } from 'lucide-react'
import { cn } from '@/lib/utils'
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from '@/components/ui/tooltip'
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select'

// Status Badge 컴포넌트
function StatusBadge({ item }: { item: any }) {
    if (item.is_manual_override) {
        return (
            <Badge className="bg-purple-500/20 text-purple-300 border-purple-500/30 text-xs">
                <UserCheck className="w-3 h-3 mr-1" />
                Manual
            </Badge>
        )
    }
    if (item.ai_refined) {
        return (
            <Badge className="bg-green-500/20 text-green-300 border-green-500/30 text-xs">
                <Sparkles className="w-3 h-3 mr-1" />
                AI Refined
            </Badge>
        )
    }
    if (!item.smiles_code) {
        return (
            <Badge className="bg-red-500/20 text-red-300 border-red-500/30 text-xs">
                <AlertTriangle className="w-3 h-3 mr-1" />
                Missing
            </Badge>
        )
    }
    return (
        <Badge variant="outline" className="text-slate-500 border-slate-600 text-xs">
            Raw
        </Badge>
    )
}

interface ReagentListTabProps {
    onSelect: (item: any) => void
    selectedId?: string
}

export function ReagentListTab({ onSelect, selectedId }: ReagentListTabProps) {
    const [data, setData] = useState<any[]>([])
    const [loading, setLoading] = useState(false)
    const [search, setSearch] = useState('')
    const [page, setPage] = useState(1)
    const [total, setTotal] = useState(0)
    const [statusFilter, setStatusFilter] = useState<string>('all')
    const limit = 20

    // Status별 카운트 (대시보드용)
    const statusCounts = {
        total: total,
        missing: data.filter(d => !d.smiles_code && !d.ai_refined).length,
        refined: data.filter(d => d.ai_refined && !d.is_manual_override).length,
        manual: data.filter(d => d.is_manual_override).length
    }

    const fetchData = async () => {
        setLoading(true)
        try {
            const API_BASE = import.meta.env.VITE_API_BASE_URL || ''
            const params = new URLSearchParams({
                page: page.toString(),
                limit: limit.toString(),
                search: search
            })
            // Add status filter
            if (statusFilter === 'missing') {
                params.append('missing_smiles', 'true')
            } else if (statusFilter === 'refined') {
                params.append('ai_refined', 'true')
            } else if (statusFilter === 'manual') {
                params.append('manual_override', 'true')
            }

            const res = await fetch(`${API_BASE}/api/library/reagents?${params}`)
            if (res.ok) {
                const json = await res.json()
                setData(json.data || [])
                setTotal(json.total || 0)
            }
        } catch (error) {
            console.error("Failed to fetch reagents", error)
        } finally {
            setLoading(false)
        }
    }

    useEffect(() => {
        const timer = setTimeout(() => {
            fetchData()
        }, 300)
        return () => clearTimeout(timer)
    }, [search, page, statusFilter])

    const totalPages = Math.ceil(total / limit)

    return (
        <div className="h-full flex flex-col bg-slate-900">
            {/* Toolbar */}
            <div className="p-4 border-b border-slate-800 space-y-3">
                {/* Status Summary Bar */}
                <div className="flex gap-2 text-xs">
                    <div className="flex items-center gap-1 px-2 py-1 bg-slate-800 rounded">
                        <span className="text-slate-400">Total:</span>
                        <span className="text-white font-medium">{total}</span>
                    </div>
                    <div className="flex items-center gap-1 px-2 py-1 bg-red-500/10 rounded border border-red-500/20">
                        <AlertTriangle className="w-3 h-3 text-red-400" />
                        <span className="text-red-300">Need Refine: {statusCounts.missing}</span>
                    </div>
                    <div className="flex items-center gap-1 px-2 py-1 bg-green-500/10 rounded border border-green-500/20">
                        <Sparkles className="w-3 h-3 text-green-400" />
                        <span className="text-green-300">AI Refined: {statusCounts.refined}</span>
                    </div>
                    <div className="flex items-center gap-1 px-2 py-1 bg-purple-500/10 rounded border border-purple-500/20">
                        <UserCheck className="w-3 h-3 text-purple-400" />
                        <span className="text-purple-300">Manual: {statusCounts.manual}</span>
                    </div>
                </div>

                {/* Search + Filter Row */}
                <div className="flex gap-2">
                    <div className="relative flex-1">
                        <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-slate-400" />
                        <Input
                            placeholder="Search reagents (Name, CAS)..."
                            value={search}
                            onChange={(e) => {
                                setSearch(e.target.value)
                                setPage(1)
                            }}
                            className="pl-9 bg-slate-800 border-slate-700 text-slate-200 focus:ring-pink-500"
                        />
                    </div>
                    <Select value={statusFilter} onValueChange={(v) => { setStatusFilter(v); setPage(1); }}>
                        <SelectTrigger className="w-[160px] bg-slate-800 border-slate-700 text-slate-200">
                            <Filter className="w-4 h-4 mr-2" />
                            <SelectValue placeholder="Filter Status" />
                        </SelectTrigger>
                        <SelectContent className="bg-slate-800 border-slate-700">
                            <SelectItem value="all">All Status</SelectItem>
                            <SelectItem value="missing">Missing SMILES</SelectItem>
                            <SelectItem value="refined">AI Refined</SelectItem>
                            <SelectItem value="manual">Manual Override</SelectItem>
                        </SelectContent>
                    </Select>
                </div>
            </div>

            {/* Table */}
            <div className="flex-1 overflow-auto">
                <Table>
                    <TableHeader className="bg-slate-900/50 sticky top-0 z-10">
                        <TableRow className="hover:bg-transparent border-slate-800">
                            <TableHead className="text-slate-400">Name</TableHead>
                            <TableHead className="text-slate-400 w-[120px]">CAS</TableHead>
                            <TableHead className="text-slate-400 w-[50px]">SMILES</TableHead>
                            <TableHead className="text-slate-400 w-[100px]">Status</TableHead>
                            <TableHead className="text-slate-400 w-[80px]">Source</TableHead>
                        </TableRow>
                    </TableHeader>
                    <TableBody>
                        {loading ? (
                            <TableRow>
                                <TableCell colSpan={5} className="h-24 text-center">
                                    <Loader2 className="w-6 h-6 animate-spin mx-auto text-pink-500" />
                                </TableCell>
                            </TableRow>
                        ) : data.length === 0 ? (
                            <TableRow>
                                <TableCell colSpan={5} className="h-24 text-center text-slate-500">
                                    No reagents found.
                                </TableCell>
                            </TableRow>
                        ) : (
                            data.map((item) => (
                                <TableRow
                                    key={item.id}
                                    className={cn(
                                        "cursor-pointer hover:bg-slate-800/50 border-slate-800 transition-colors",
                                        selectedId === item.id && "bg-pink-900/20 hover:bg-pink-900/30 border-l-2 border-l-pink-500",
                                        // SMILES 없고 정제 안 된 항목은 행 배경 강조
                                        !item.smiles_code && !item.ai_refined && "bg-red-950/20"
                                    )}
                                    onClick={() => onSelect(item)}
                                >
                                    <TableCell className="font-medium text-slate-200">
                                        <div className="truncate max-w-[200px]" title={item.product_name || item.name}>
                                            {item.product_name || item.name || '-'}
                                        </div>
                                    </TableCell>
                                    <TableCell className="text-slate-400 font-mono text-xs">
                                        {item.cas_number || '-'}
                                    </TableCell>
                                    <TableCell>
                                        {item.smiles_code ? (
                                            <TooltipProvider>
                                                <Tooltip>
                                                    <TooltipTrigger>
                                                        <CheckCircle2 className="w-4 h-4 text-green-500" />
                                                    </TooltipTrigger>
                                                    <TooltipContent>
                                                        <p className="font-mono text-xs max-w-xs break-all">{item.smiles_code}</p>
                                                    </TooltipContent>
                                                </Tooltip>
                                            </TooltipProvider>
                                        ) : (
                                            <TooltipProvider>
                                                <Tooltip>
                                                    <TooltipTrigger>
                                                        <AlertTriangle className="w-4 h-4 text-red-500" />
                                                    </TooltipTrigger>
                                                    <TooltipContent>
                                                        <p>Missing SMILES - Needs refinement</p>
                                                    </TooltipContent>
                                                </Tooltip>
                                            </TooltipProvider>
                                        )}
                                    </TableCell>
                                    <TableCell>
                                        <StatusBadge item={item} />
                                    </TableCell>
                                    <TableCell className="text-slate-400 text-xs">
                                        {item.source_name || 'Unknown'}
                                    </TableCell>
                                </TableRow>
                            ))
                        )}
                    </TableBody>
                </Table>
            </div>

            {/* Pagination */}
            <div className="p-2 border-t border-slate-800 flex items-center justify-between text-xs text-slate-400">
                <span>Total {total} items</span>
                <div className="flex gap-1">
                    <Button
                        variant="ghost"
                        size="icon"
                        disabled={page <= 1}
                        onClick={() => setPage(p => p - 1)}
                        className="h-8 w-8 hover:bg-slate-800 hover:text-white"
                    >
                        <ChevronLeft className="w-4 h-4" />
                    </Button>
                    <span className="flex items-center px-2">
                        Page {page} of {totalPages || 1}
                    </span>
                    <Button
                        variant="ghost"
                        size="icon"
                        disabled={page >= totalPages}
                        onClick={() => setPage(p => p + 1)}
                        className="h-8 w-8 hover:bg-slate-800 hover:text-white"
                    >
                        <ChevronRight className="w-4 h-4" />
                    </Button>
                </div>
            </div>
        </div>
    )
}
