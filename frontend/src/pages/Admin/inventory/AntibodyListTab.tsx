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
import { Search, Loader2, ChevronLeft, ChevronRight, FlaskConical } from 'lucide-react'
import { cn } from '@/lib/utils'

interface AntibodyListTabProps {
    onSelect: (item: any) => void
    selectedId?: string
}

export function AntibodyListTab({ onSelect, selectedId }: AntibodyListTabProps) {
    const [data, setData] = useState<any[]>([])
    const [loading, setLoading] = useState(false)
    const [search, setSearch] = useState('')
    const [page, setPage] = useState(1)
    const [total, setTotal] = useState(0)
    const limit = 20

    const fetchData = async () => {
        setLoading(true)
        try {
            const API_BASE = import.meta.env.VITE_API_BASE_URL || ''
            const params = new URLSearchParams({
                page: page.toString(),
                limit: limit.toString(),
                search: search
            })
            const res = await fetch(`${API_BASE}/api/library/antibodies?${params}`)
            if (res.ok) {
                const json = await res.json()
                setData(json.data || [])
                setTotal(json.total || 0)
            }
        } catch (error) {
            console.error("Failed to fetch antibodies", error)
        } finally {
            setLoading(false)
        }
    }

    useEffect(() => {
        const timer = setTimeout(() => {
            fetchData()
        }, 300)
        return () => clearTimeout(timer)
    }, [search, page])

    const totalPages = Math.ceil(total / limit)

    return (
        <div className="h-full flex flex-col bg-slate-900">
            {/* Toolbar */}
            <div className="p-4 border-b border-slate-800 flex gap-2">
                <div className="relative flex-1">
                    <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-slate-400" />
                    <Input
                        placeholder="Search antibodies..."
                        value={search}
                        onChange={(e) => {
                            setSearch(e.target.value)
                            setPage(1)
                        }}
                        className="pl-9 bg-slate-800 border-slate-700 text-slate-200 focus:ring-blue-500"
                    />
                </div>
            </div>

            {/* Table */}
            <div className="flex-1 overflow-auto">
                <Table>
                    <TableHeader className="bg-slate-900/50 sticky top-0 z-10">
                        <TableRow className="hover:bg-transparent border-slate-800">
                            <TableHead className="text-slate-400 w-[100px]">Cat No</TableHead>
                            <TableHead className="text-slate-400">Name</TableHead>
                            <TableHead className="text-slate-400">Target</TableHead>
                            <TableHead className="text-slate-400">Host</TableHead>
                            <TableHead className="text-slate-400 w-[80px]">Source</TableHead>
                        </TableRow>
                    </TableHeader>
                    <TableBody>
                        {loading ? (
                            <TableRow>
                                <TableCell colSpan={5} className="h-24 text-center">
                                    <Loader2 className="w-6 h-6 animate-spin mx-auto text-blue-500" />
                                </TableCell>
                            </TableRow>
                        ) : data.length === 0 ? (
                            <TableRow>
                                <TableCell colSpan={5} className="h-24 text-center text-slate-500">
                                    No antibodies found.
                                </TableCell>
                            </TableRow>
                        ) : (
                            data.map((item) => (
                                <TableRow
                                    key={item.id}
                                    className={cn(
                                        "cursor-pointer hover:bg-slate-800/50 border-slate-800 transition-colors",
                                        selectedId === item.id && "bg-blue-900/20 hover:bg-blue-900/30 border-l-2 border-l-blue-500"
                                    )}
                                    onClick={() => onSelect(item)}
                                >
                                    <TableCell className="font-mono text-xs text-blue-400">
                                        {item.cat_no || '-'}
                                    </TableCell>
                                    <TableCell className="font-medium text-slate-200">
                                        <div className="truncate max-w-[200px]" title={item.product_name}>
                                            {item.product_name || item.name || '-'}
                                        </div>
                                    </TableCell>
                                    <TableCell className="text-slate-300">
                                        {item.related_disease || item.target || '-'}
                                    </TableCell>
                                    <TableCell className="text-slate-400">
                                        {item.host_species || '-'}
                                    </TableCell>
                                    <TableCell>
                                        <Badge variant="outline" className="text-xs border-slate-600 text-slate-400">
                                            {item.source_name || 'Unknown'}
                                        </Badge>
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
