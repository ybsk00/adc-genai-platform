import { useState, useEffect } from 'react'
import {
    Table,
    TableBody,
    TableCell,
    TableHead,
    TableHeader,
    TableRow,
} from "@/components/ui/table"
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Badge } from '@/components/ui/badge'
import {
    Search,
    RefreshCw,
    Edit2,
    Loader2,
    ChevronLeft,
    ChevronRight,
    Check,
    FileText
} from 'lucide-react'
import { API_BASE_URL } from '@/lib/api'
import { toast } from 'sonner'
import {
    Dialog,
    DialogContent,
    DialogHeader,
    DialogTitle,
    DialogFooter,
} from "@/components/ui/dialog"
import { Label } from "@/components/ui/label"
import { Textarea } from "@/components/ui/textarea"

interface KnowledgeItem {
    id: string
    title: string
    content: string | null
    summary: string | null
    relevance_score: number | null
    rag_status: string
    created_at: string
}

interface KnowledgeTableProps {
    missingDataOnly: boolean
}

export function KnowledgeTable({ missingDataOnly }: KnowledgeTableProps) {
    const [data, setData] = useState<KnowledgeItem[]>([])
    const [loading, setLoading] = useState(true)
    const [search, setSearch] = useState('')
    const [page, setPage] = useState(0)
    const [total, setTotal] = useState(0)
    const limit = 10

    // Edit State
    const [editingItem, setEditingItem] = useState<KnowledgeItem | null>(null)
    const [editForm, setEditForm] = useState<Partial<KnowledgeItem>>({})
    const [saving, setSaving] = useState(false)

    const fetchData = async () => {
        setLoading(true)
        try {
            const params = new URLSearchParams({
                limit: limit.toString(),
                offset: (page * limit).toString(),
                search: search,
                missing_data_only: missingDataOnly ? 'true' : 'false'
            })
            
            const res = await fetch(`${API_BASE_URL}/api/admin/inventory/knowledge?${params}`)
            if (res.ok) {
                const json = await res.json()
                setData(json.data)
                setTotal(json.total)
            }
        } catch (error) {
            toast.error("Failed to fetch knowledge base")
        } finally {
            setLoading(false)
        }
    }

    useEffect(() => {
        fetchData()
    }, [page, search, missingDataOnly])

    const openEdit = (item: KnowledgeItem) => {
        setEditingItem(item)
        setEditForm({
            title: item.title,
            summary: item.summary,
            content: item.content
        })
    }

    const handleSave = async () => {
        if (!editingItem) return
        setSaving(true)
        try {
            const res = await fetch(`${API_BASE_URL}/api/admin/inventory/knowledge_base/${editingItem.id}`, {
                method: 'PATCH',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ updates: editForm })
            })
            
            if (res.ok) {
                toast.success("Item updated")
                setEditingItem(null)
                fetchData() // Refresh
            } else {
                toast.error("Failed to update")
            }
        } catch (e) {
            toast.error("Error saving")
        } finally {
            setSaving(false)
        }
    }

    return (
        <div className="space-y-4">
            <div className="flex gap-2 justify-between items-center mb-4">
                <div className="relative w-72">
                    <Search className="absolute left-3 top-3 w-4 h-4 text-slate-400" />
                    <Input 
                        placeholder="Search articles..." 
                        className="pl-9 bg-slate-950 border-slate-700"
                        value={search}
                        onChange={(e) => setSearch(e.target.value)}
                    />
                </div>
                <div className="flex items-center gap-2">
                    <span className="text-xs text-slate-500">
                        Total {total.toLocaleString()} items
                    </span>
                    <Button variant="outline" size="sm" onClick={fetchData} className="border-slate-700">
                        <RefreshCw className="w-3 h-3 mr-2" /> Refresh
                    </Button>
                </div>
            </div>
            
            <div className="rounded-md border border-slate-800 overflow-hidden">
                <Table>
                    <TableHeader className="bg-slate-950">
                        <TableRow>
                            <TableHead className="w-[400px]">Title</TableHead>
                            <TableHead>Summary</TableHead>
                            <TableHead>Score</TableHead>
                            <TableHead>RAG Status</TableHead>
                            <TableHead className="text-right">Actions</TableHead>
                        </TableRow>
                    </TableHeader>
                    <TableBody>
                        {loading ? (
                            <TableRow>
                                <TableCell colSpan={5} className="h-24 text-center">
                                    <Loader2 className="w-6 h-6 animate-spin mx-auto text-slate-500" />
                                </TableCell>
                            </TableRow>
                        ) : data.length === 0 ? (
                            <TableRow>
                                <TableCell colSpan={5} className="h-24 text-center text-slate-500">
                                    No items found.
                                </TableCell>
                            </TableRow>
                        ) : (
                            data.map((item) => (
                                <TableRow key={item.id} className="hover:bg-slate-800/50 transition-colors group border-slate-800/50">
                                    <TableCell className="font-medium text-slate-200 align-top py-4">
                                        <div className="flex gap-2">
                                            <FileText className="w-4 h-4 text-slate-500 mt-1 shrink-0" />
                                            <div>
                                                {item.title}
                                                <div className="text-[10px] text-slate-500 font-mono mt-1">{item.created_at.slice(0, 10)}</div>
                                            </div>
                                        </div>
                                    </TableCell>
                                    <TableCell className="text-sm text-slate-400 align-top py-4 max-w-[300px]">
                                        <div className="line-clamp-2">
                                            {item.summary || <span className="text-slate-600 italic">No summary generated</span>}
                                        </div>
                                    </TableCell>
                                    <TableCell className="align-top py-4">
                                        {item.relevance_score !== null ? (
                                            <Badge variant="outline" className={`
                                                ${(item.relevance_score || 0) >= 0.7 ? 'border-green-500/50 text-green-400' : 'border-slate-700 text-slate-500'}
                                            `}>
                                                {item.relevance_score.toFixed(2)}
                                            </Badge>
                                        ) : (
                                            <span className="text-slate-600 text-xs">-</span>
                                        )}
                                    </TableCell>
                                    <TableCell className="align-top py-4">
                                        {item.rag_status === 'processed' ? (
                                            <Badge className="bg-blue-500/10 text-blue-400 hover:bg-blue-500/20 text-[10px] px-1 py-0 border-0">
                                                Indexed
                                            </Badge>
                                        ) : (
                                            <Badge variant="outline" className="text-slate-500 border-slate-700 text-[10px]">
                                                {item.rag_status}
                                            </Badge>
                                        )}
                                    </TableCell>
                                    <TableCell className="text-right align-top py-4">
                                        <Button 
                                            size="icon" variant="ghost" 
                                            className="h-7 w-7 text-blue-400 hover:bg-blue-400/10 hover:text-blue-300 opacity-60 group-hover:opacity-100" 
                                            title="Edit"
                                            onClick={() => openEdit(item)}
                                        >
                                            <Edit2 className="w-3 h-3" />
                                        </Button>
                                    </TableCell>
                                </TableRow>
                            ))
                        )}
                    </TableBody>
                </Table>
            </div>

            {/* Pagination */}
            <div className="flex justify-end items-center gap-2 mt-4">
                <Button variant="outline" size="sm" onClick={() => setPage(p => Math.max(0, p - 1))} disabled={page === 0}>
                    <ChevronLeft className="w-4 h-4" />
                </Button>
                <span className="text-xs text-slate-400">
                    Page {page + 1} of {Math.ceil(total / limit)}
                </span>
                <Button variant="outline" size="sm" onClick={() => setPage(p => p + 1)} disabled={(page + 1) * limit >= total}>
                    <ChevronRight className="w-4 h-4" />
                </Button>
            </div>

            {/* Edit Dialog */}
            <Dialog open={!!editingItem} onOpenChange={(open) => !open && setEditingItem(null)}>
                <DialogContent className="bg-slate-900 border-slate-700 text-slate-100 sm:max-w-[600px]">
                    <DialogHeader>
                        <DialogTitle>Edit Knowledge Item</DialogTitle>
                    </DialogHeader>
                    <div className="grid gap-4 py-4">
                        <div className="grid grid-cols-4 items-center gap-4">
                            <Label htmlFor="title" className="text-right text-slate-400">Title</Label>
                            <Textarea 
                                id="title" 
                                value={editForm.title || ''} 
                                onChange={(e) => setEditForm({...editForm, title: e.target.value})}
                                className="col-span-3 bg-slate-950 border-slate-700" 
                            />
                        </div>
                        <div className="grid grid-cols-4 items-center gap-4">
                            <Label htmlFor="summary" className="text-right text-slate-400">Summary</Label>
                            <Textarea 
                                id="summary" 
                                value={editForm.summary || ''} 
                                onChange={(e) => setEditForm({...editForm, summary: e.target.value})}
                                className="col-span-3 bg-slate-950 border-slate-700 min-h-[100px]" 
                            />
                        </div>
                        <div className="grid grid-cols-4 items-center gap-4">
                            <Label htmlFor="content" className="text-right text-slate-400">Content</Label>
                            <Textarea 
                                id="content" 
                                value={editForm.content || ''} 
                                onChange={(e) => setEditForm({...editForm, content: e.target.value})}
                                className="col-span-3 bg-slate-950 border-slate-700 min-h-[150px]" 
                            />
                        </div>
                    </div>
                    <DialogFooter>
                        <Button variant="ghost" onClick={() => setEditingItem(null)}>Cancel</Button>
                        <Button onClick={handleSave} disabled={saving} className="bg-blue-600 hover:bg-blue-700 text-white">
                            {saving ? <Loader2 className="w-4 h-4 animate-spin mr-2" /> : <Check className="w-4 h-4 mr-2" />}
                            Save & Re-Index
                        </Button>
                    </DialogFooter>
                </DialogContent>
            </Dialog>
        </div>
    )
}
