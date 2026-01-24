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
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import {
    Search,
    RefreshCw,
    Sparkles,
    CheckCircle2,
    XCircle,
    Edit2,
    Save,
    Loader2
} from 'lucide-react'
import { API_BASE_URL } from '@/lib/api'
import { toast } from 'sonner'
import { Card, CardContent } from '@/components/ui/card'
import {
    Dialog,
    DialogContent,
    DialogHeader,
    DialogTitle,
    DialogFooter
} from '@/components/ui/dialog'
import { Label } from "@/components/ui/label"

interface GoldenSetItem {
    id: string
    drug_name: string
    name?: string // Add name as optional
    target: string | null
    outcome_type: string | null
    failure_reason: string | null
    created_at: string
    is_ai_extracted: boolean
    properties: any
}

function PromotedTable({ tab }: { tab: 'success' | 'failure' }) {
    const [data, setData] = useState<GoldenSetItem[]>([])
    const [loading, setLoading] = useState(true)
    const [search, setSearch] = useState('')
    
    // Edit
    const [editingItem, setEditingItem] = useState<GoldenSetItem | null>(null)
    const [editForm, setEditForm] = useState<Partial<GoldenSetItem>>({})
    const [saving, setSaving] = useState(false)

    const fetchData = async () => {
        setLoading(true)
        try {
            const res = await fetch(`${API_BASE_URL}/api/admin/goldenset/library?limit=50&tab=${tab}&search=${search}`)
            if (res.ok) {
                const json = await res.json()
                setData(json.data)
            }
        } catch (error) {
            toast.error("Failed to fetch library")
        } finally {
            setLoading(false)
        }
    }

    useEffect(() => {
        fetchData()
    }, [tab, search])

    const handleSave = async () => {
        if (!editingItem) return
        setSaving(true)
        try {
            // Update table='golden_set_library'
            const res = await fetch(`${API_BASE_URL}/api/admin/inventory/golden_set_library/${editingItem.id}`, {
                method: 'PATCH',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ 
                    updates: { 
                        name: editForm.drug_name,
                        outcome_type: editForm.outcome_type,
                        failure_reason: editForm.failure_reason,
                        properties: { ...editingItem.properties, ...editForm.properties }
                    } 
                })
            })
            
            if (res.ok) {
                toast.success("Updated successfully")
                if (editForm.outcome_type === 'Success' && tab === 'failure') {
                    toast.success("Item moved to Success tab!")
                }
                setEditingItem(null)
                fetchData()
            } else {
                toast.error("Update failed")
            }
        } catch (e) {
            toast.error("Error saving")
        } finally {
            setSaving(false)
        }
    }

    return (
        <div className="space-y-4">
            <div className="flex justify-between items-center mb-4">
                <div className="relative w-72">
                    <Search className="absolute left-3 top-3 w-4 h-4 text-slate-400" />
                    <Input 
                        placeholder="Search..." 
                        className="pl-9 bg-slate-950 border-slate-700"
                        value={search}
                        onChange={(e) => setSearch(e.target.value)}
                    />
                </div>
                <Button variant="outline" size="sm" onClick={fetchData} className="border-slate-700">
                    <RefreshCw className="w-3 h-3 mr-2" /> Refresh
                </Button>
            </div>

            <div className="rounded-md border border-slate-800 overflow-hidden">
                <Table>
                    <TableHeader className="bg-slate-950">
                        <TableRow>
                            <TableHead className="w-[300px]">Drug Name</TableHead>
                            <TableHead>Target</TableHead>
                            <TableHead>Outcome</TableHead>
                            <TableHead>Reason</TableHead>
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
                                <TableRow key={item.id} className="hover:bg-slate-800/50 border-slate-800/50">
                                    <TableCell className="font-medium text-slate-200">
                                        {item.drug_name || item['name']}
                                    </TableCell>
                                    <TableCell>{item.target || item.properties?.target}</TableCell>
                                    <TableCell>
                                        <Badge variant="outline" className={
                                            item.outcome_type === 'Success' 
                                            ? "border-green-500/50 text-green-400" 
                                            : "border-red-500/50 text-red-400"
                                        }>
                                            {item.outcome_type || 'Unknown'}
                                        </Badge>
                                    </TableCell>
                                    <TableCell className="text-xs text-slate-400 max-w-[200px] truncate">
                                        {item.failure_reason}
                                    </TableCell>
                                    <TableCell className="text-right">
                                        <Button 
                                            size="icon" variant="ghost" 
                                            className="h-7 w-7 text-blue-400 hover:bg-blue-400/10"
                                            onClick={() => {
                                                setEditingItem(item)
                                                setEditForm({
                                                    drug_name: item.drug_name || item['name'],
                                                    outcome_type: item.outcome_type,
                                                    failure_reason: item.failure_reason
                                                })
                                            }}
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

             {/* Edit Dialog */}
             <Dialog open={!!editingItem} onOpenChange={(open) => !open && setEditingItem(null)}>
                <DialogContent className="bg-slate-900 border-slate-700 text-slate-100">
                    <DialogHeader>
                        <DialogTitle>Edit Golden Set Item</DialogTitle>
                    </DialogHeader>
                    <div className="grid gap-4 py-4">
                        <div className="grid grid-cols-4 items-center gap-4">
                            <Label htmlFor="name" className="text-right text-slate-400">Name</Label>
                            <Input 
                                id="name" 
                                value={editForm.drug_name || ''} 
                                onChange={(e) => setEditForm({...editForm, drug_name: e.target.value})}
                                className="col-span-3 bg-slate-950 border-slate-700" 
                            />
                        </div>
                        <div className="grid grid-cols-4 items-center gap-4">
                            <Label htmlFor="outcome" className="text-right text-slate-400">Outcome</Label>
                            <div className="col-span-3 flex gap-2">
                                <Button 
                                    size="sm" 
                                    variant={editForm.outcome_type === 'Success' ? 'default' : 'outline'}
                                    className={editForm.outcome_type === 'Success' ? 'bg-green-600' : 'border-slate-700'}
                                    onClick={() => setEditForm({...editForm, outcome_type: 'Success', failure_reason: null})}
                                >
                                    Success
                                </Button>
                                <Button 
                                    size="sm" 
                                    variant={editForm.outcome_type === 'Failure' ? 'destructive' : 'outline'}
                                    onClick={() => setEditForm({...editForm, outcome_type: 'Failure'})}
                                >
                                    Failure
                                </Button>
                            </div>
                        </div>
                        {editForm.outcome_type !== 'Success' && (
                            <div className="grid grid-cols-4 items-center gap-4">
                                <Label htmlFor="reason" className="text-right text-slate-400">Reason</Label>
                                <Input 
                                    id="reason" 
                                    value={editForm.failure_reason || ''} 
                                    onChange={(e) => setEditForm({...editForm, failure_reason: e.target.value})}
                                    className="col-span-3 bg-slate-950 border-slate-700" 
                                />
                            </div>
                        )}
                    </div>
                    <DialogFooter>
                        <Button variant="ghost" onClick={() => setEditingItem(null)}>Cancel</Button>
                        <Button onClick={handleSave} disabled={saving} className="bg-blue-600 hover:bg-blue-700 text-white">
                            {saving ? <Loader2 className="w-4 h-4 animate-spin mr-2" /> : <Save className="w-4 h-4 mr-2" />}
                            Save Changes
                        </Button>
                    </DialogFooter>
                </DialogContent>
            </Dialog>
        </div>
    )
}

export function GoldenSetLibraryPage() {
    return (
        <div className="space-y-6">
            <div>
                <h1 className="text-3xl font-bold text-white flex items-center gap-3">
                    <Sparkles className="w-8 h-8 text-yellow-500" />
                    Golden Set Library
                </h1>
                <p className="text-slate-400 mt-2">
                    Repository of all approved Golden Set items.
                </p>
            </div>
            
            <Card className="bg-slate-900 border-slate-800">
                <CardContent className="p-6">
                    <Tabs defaultValue="success" className="w-full">
                        <TabsList className="grid w-full grid-cols-2 bg-slate-800 mb-6">
                            <TabsTrigger value="success" className="data-[state=active]:bg-green-900/50 data-[state=active]:text-green-200">
                                <CheckCircle2 className="w-4 h-4 mr-2" />
                                Success (Approved)
                            </TabsTrigger>
                            <TabsTrigger value="failure" className="data-[state=active]:bg-red-900/50 data-[state=active]:text-red-200">
                                <XCircle className="w-4 h-4 mr-2" />
                                Failure / Pending (Needs Review)
                            </TabsTrigger>
                        </TabsList>

                        <TabsContent value="success">
                            <PromotedTable tab="success" />
                        </TabsContent>
                        <TabsContent value="failure">
                            <PromotedTable tab="failure" />
                        </TabsContent>
                    </Tabs>
                </CardContent>
            </Card>
        </div>
    )
}
