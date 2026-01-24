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
    ArrowUpRight,
    Edit2,
    Sparkles,
    Loader2,
    ChevronLeft,
    ChevronRight,
    Check
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

interface CommercialItem {
    id: string
    product_name: string
    target: string | null
    category: string | null
    smiles_code: string | null
    summary: string | null
    ai_refined: boolean
    source_name: string
    crawled_at: string
}

interface CommercialTableProps {
    missingDataOnly: boolean
}

export function CommercialTable({ missingDataOnly }: CommercialTableProps) {
    const [data, setData] = useState<CommercialItem[]>([])
    const [loading, setLoading] = useState(true)
    const [search, setSearch] = useState('')
    const [page, setPage] = useState(0)
    const [total, setTotal] = useState(0)
    const limit = 10

    // Edit State
    const [editingItem, setEditingItem] = useState<CommercialItem | null>(null)
    const [editForm, setEditForm] = useState<Partial<CommercialItem>>({})
    const [saving, setSaving] = useState(false)

    // Promote State
    const [promoting, setPromoting] = useState<string | null>(null)

    const fetchData = async () => {
        setLoading(true)
        try {
            const params = new URLSearchParams({
                limit: limit.toString(),
                offset: (page * limit).toString(),
                search: search,
                missing_data_only: missingDataOnly ? 'true' : 'false'
            })
            
            const res = await fetch(`${API_BASE_URL}/api/admin/inventory/commercial?${params}`)
            if (res.ok) {
                const json = await res.json()
                setData(json.data)
                setTotal(json.total)
            }
        } catch (error) {
            toast.error("Failed to fetch data")
        } finally {
            setLoading(false)
        }
    }

    useEffect(() => {
        fetchData()
    }, [page, search, missingDataOnly])

    const handlePromote = async (id: string) => {
        if (!confirm("Promote this item to Golden Set?")) return
        
        setPromoting(id)
        try {
            const res = await fetch(`${API_BASE_URL}/api/admin/inventory/promote/${id}`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ target_table: 'golden_set_library' })
            })
            
            if (res.ok) {
                toast.success("Successfully promoted to Golden Set!")
            } else {
                toast.error("Failed to promote")
            }
        } catch (e) {
            toast.error("Error promoting item")
        } finally {
            setPromoting(null)
        }
    }

    const openEdit = (item: CommercialItem) => {
        setEditingItem(item)
        setEditForm({
            product_name: item.product_name,
            target: item.target,
            smiles_code: item.smiles_code,
            summary: item.summary
        })
    }

    const handleSave = async () => {
        if (!editingItem) return
        setSaving(true)
        try {
            const res = await fetch(`${API_BASE_URL}/api/admin/inventory/commercial_reagents/${editingItem.id}`, {
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
                        placeholder="Search reagents..." 
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
                            <TableHead className="w-[300px]">Product Name</TableHead>
                            <TableHead>Target</TableHead>
                            <TableHead>Source</TableHead>
                            <TableHead>Structure (SMILES)</TableHead>
                            <TableHead>AI Status</TableHead>
                            <TableHead className="text-right">Actions</TableHead>
                        </TableRow>
                    </TableHeader>
                    <TableBody>
                        {loading ? (
                            <TableRow>
                                <TableCell colSpan={6} className="h-24 text-center">
                                    <Loader2 className="w-6 h-6 animate-spin mx-auto text-slate-500" />
                                </TableCell>
                            </TableRow>
                        ) : data.length === 0 ? (
                            <TableRow>
                                <TableCell colSpan={6} className="h-24 text-center text-slate-500">
                                    No items found.
                                </TableCell>
                            </TableRow>
                        ) : (
                            data.map((item) => (
                                <TableRow key={item.id} className="hover:bg-slate-800/50 transition-colors group border-slate-800/50">
                                    <TableCell className="font-medium text-slate-200">
                                        {item.product_name}
                                        <div className="text-[10px] text-slate-500 font-mono mt-0.5">{item.id.slice(0, 8)}...</div>
                                    </TableCell>
                                    <TableCell>
                                        {item.target ? (
                                            <Badge variant="outline" className="bg-slate-800 text-blue-300 border-slate-600">
                                                {item.target}
                                            </Badge>
                                        ) : (
                                            <span className="text-slate-600 text-xs italic">Unknown</span>
                                        )}
                                    </TableCell>
                                    <TableCell className="text-xs text-slate-400">
                                        {item.source_name}
                                    </TableCell>
                                    <TableCell>
                                        {item.smiles_code ? (
                                            <div className="text-xs font-mono text-slate-400 truncate max-w-[120px]" title={item.smiles_code}>
                                                {item.smiles_code}
                                            </div>
                                        ) : (
                                            <span className="text-slate-600 text-xs italic">Missing</span>
                                        )}
                                    </TableCell>
                                    <TableCell>
                                        {item.ai_refined ? (
                                            <Badge className="bg-green-500/10 text-green-500 hover:bg-green-500/20 text-[10px] px-1 py-0 border-0">
                                                Refined
                                            </Badge>
                                        ) : (
                                            <Badge variant="outline" className="text-slate-500 border-slate-700 text-[10px]">
                                                Pending
                                            </Badge>
                                        )}
                                    </TableCell>
                                    <TableCell className="text-right">
                                        <div className="flex justify-end gap-1 opacity-60 group-hover:opacity-100 transition-opacity">
                                            <Button 
                                                size="icon" variant="ghost" 
                                                className="h-7 w-7 text-yellow-500 hover:bg-yellow-500/10 hover:text-yellow-400" 
                                                title="Promote to Golden Set"
                                                onClick={() => handlePromote(item.id)}
                                                disabled={!!promoting}
                                            >
                                                {promoting === item.id ? <Loader2 className="w-3 h-3 animate-spin"/> : <ArrowUpRight className="w-3 h-3" />}
                                            </Button>
                                            <Button 
                                                size="icon" variant="ghost" 
                                                className="h-7 w-7 text-blue-400 hover:bg-blue-400/10 hover:text-blue-300" 
                                                title="Edit"
                                                onClick={() => openEdit(item)}
                                            >
                                                <Edit2 className="w-3 h-3" />
                                            </Button>
                                        </div>
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
                <DialogContent className="bg-slate-900 border-slate-700 text-slate-100 sm:max-w-[500px]">
                    <DialogHeader>
                        <DialogTitle>Edit Reagent</DialogTitle>
                    </DialogHeader>
                    <div className="grid gap-4 py-4">
                        <div className="grid grid-cols-4 items-center gap-4">
                            <Label htmlFor="name" className="text-right text-slate-400">Name</Label>
                            <Input 
                                id="name" 
                                value={editForm.product_name || ''} 
                                onChange={(e) => setEditForm({...editForm, product_name: e.target.value})}
                                className="col-span-3 bg-slate-950 border-slate-700" 
                            />
                        </div>
                        <div className="grid grid-cols-4 items-center gap-4">
                            <Label htmlFor="target" className="text-right text-slate-400">Target</Label>
                            <Input 
                                id="target" 
                                value={editForm.target || ''} 
                                onChange={(e) => setEditForm({...editForm, target: e.target.value})}
                                className="col-span-3 bg-slate-950 border-slate-700" 
                            />
                        </div>
                        <div className="grid grid-cols-4 items-center gap-4">
                            <Label htmlFor="smiles" className="text-right text-slate-400">SMILES</Label>
                            <Input 
                                id="smiles" 
                                value={editForm.smiles_code || ''} 
                                onChange={(e) => setEditForm({...editForm, smiles_code: e.target.value})}
                                className="col-span-3 bg-slate-950 border-slate-700 font-mono text-xs" 
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
                    </div>
                    <DialogFooter>
                        <Button variant="ghost" onClick={() => setEditingItem(null)}>Cancel</Button>
                        <Button onClick={handleSave} disabled={saving} className="bg-blue-600 hover:bg-blue-700 text-white">
                            {saving ? <Loader2 className="w-4 h-4 animate-spin mr-2" /> : <Check className="w-4 h-4 mr-2" />}
                            Save Changes
                        </Button>
                    </DialogFooter>
                </DialogContent>
            </Dialog>
        </div>
    )
}
