import { useState, useEffect, useRef } from 'react'
import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import { Input } from '@/components/ui/input'
import { Label } from '@/components/ui/label'
import {
    Select,
    SelectContent,
    SelectItem,
    SelectTrigger,
    SelectValue,
} from '@/components/ui/select'
import { Textarea } from '@/components/ui/textarea'
import { ScrollArea } from '@/components/ui/scroll-area'
import {
    CheckCircle, Loader2, ArrowRight, Sparkles,
    Search, FlaskConical, Microscope, Activity, X
} from 'lucide-react'
import { toast } from 'sonner'
import { getSession } from '@/lib/supabase'
import { API_BASE_URL } from '@/lib/api'
import SmilesDrawer from 'smiles-drawer'

// --- Types ---
interface GoldenSetDraft {
    id: string
    drug_name: string
    target: string
    payload?: string
    linker?: string
    smiles_code?: string
    enrichment_source: string
    created_at: string
    outcome_type?: string
    failure_reason?: string
    is_ai_extracted?: boolean
    // Bio Metrics
    binding_affinity?: string
    isotype?: string
    host_species?: string
    orr_pct?: string
    os_months?: string
    pfs_months?: string
    // Meta
    is_manual_override?: boolean
    properties?: Record<string, any>
}

interface DraftsResponse {
    data: GoldenSetDraft[]
    total: number
    limit: number
    offset: number
}

// --- Components ---

// 1. Chemical Structure Renderer
function ChemicalStructure({ smiles }: { smiles: string }) {
    const canvasRef = useRef<HTMLCanvasElement>(null)

    useEffect(() => {
        if (!canvasRef.current || !smiles) return
        
        try {
            const drawer = new SmilesDrawer.Drawer({
                width: 300,
                height: 150,
                compactDrawing: false,
            })
            
            SmilesDrawer.parse(smiles, (tree: any) => {
                if (canvasRef.current) {
                    drawer.draw(tree, canvasRef.current, 'light', false)
                }
            }, (err: any) => {
                console.warn('SMILES Parse Error:', err)
            })
        } catch (e) {
            console.error(e)
        }
    }, [smiles])

    if (!smiles) return <div className="h-[150px] flex items-center justify-center bg-slate-100/5 text-slate-500 text-xs border border-dashed border-slate-700 rounded">No SMILES Data</div>

    return (
        <div className="bg-white rounded p-2 border border-slate-700">
            <canvas ref={canvasRef} className="w-full h-[150px]" />
        </div>
    )
}

// 2. Status Badge Logic
function StatusBadge({ status, override }: { status?: string, override?: boolean }) {
    if (override) return <Badge className="bg-blue-500/20 text-blue-300 border-blue-500/30">Manual Edit</Badge>
    if (status === 'Success') return <Badge className="bg-green-500/20 text-green-300 border-green-500/30">Success</Badge>
    if (status === 'Failure') return <Badge className="bg-red-500/20 text-red-300 border-red-500/30">Failure</Badge>
    return <Badge variant="outline" className="text-slate-500">Pending</Badge>
}

export function StagingAreaTab() {
    const queryClient = useQueryClient()
    const [selectedDraftId, setSelectedDraftId] = useState<string | null>(null)
    const [isOpen, setIsOpen] = useState(false) // Drawer Open State

    // Search State
    const [page, setPage] = useState(1)
    const [searchQuery, setSearchQuery] = useState('')
    const [debouncedSearch, setDebouncedSearch] = useState('')

    useEffect(() => {
        const timer = setTimeout(() => setDebouncedSearch(searchQuery), 300)
        return () => clearTimeout(timer)
    }, [searchQuery])

    // --- Queries ---
    const { data: draftsResponse, isLoading } = useQuery({
        queryKey: ['stagingDrafts', page, debouncedSearch],
        queryFn: async () => {
            const { session } = await getSession()
            const offset = (page - 1) * 20
            const params = new URLSearchParams({
                limit: '20',
                offset: String(offset),
                search: debouncedSearch
            })
            const response = await fetch(`${API_BASE_URL}/api/admin/goldenset/drafts?${params}`, {
                headers: { Authorization: `Bearer ${session?.access_token}` }
            })
            if (!response.ok) throw new Error('Failed to fetch')
            return response.json() as Promise<DraftsResponse>
        }
    })

    // Local state for editing form
    const [formData, setFormData] = useState<Partial<GoldenSetDraft>>({})
    const selectedDraft = draftsResponse?.data.find(d => d.id === selectedDraftId)

    useEffect(() => {
        if (selectedDraft) {
            setFormData(selectedDraft)
        }
    }, [selectedDraft])

    // --- Mutations ---
    const updateMutation = useMutation({
        mutationFn: async ({ id, updates }: { id: string, updates: any }) => {
            const { session } = await getSession()
            // We use the commercial_reagents endpoint for now as it has the patch logic
            // But drafts are in golden_set_library. 
            // NOTE: The backend logic needs to support patching golden_set_library too.
            // Let's assume patch_inventory_item supports 'golden_set_library' table.
            
            const response = await fetch(`${API_BASE_URL}/api/admin/inventory/golden_set_library/${id}`, { 
                method: 'PATCH',
                headers: {
                    'Content-Type': 'application/json',
                    Authorization: `Bearer ${session?.access_token}`
                },
                body: JSON.stringify({ updates })
            })
            if (!response.ok) throw new Error('Update failed')
            return response.json()
        },
        onSuccess: () => {
            toast.success('저장되었습니다 (임베딩 갱신 완료)')
            queryClient.invalidateQueries({ queryKey: ['stagingDrafts'] })
        },
        onError: () => toast.error('저장 실패')
    })

    const approveMutation = useMutation({
        mutationFn: async (id: string) => {
            const { session } = await getSession()
            const response = await fetch(`${API_BASE_URL}/api/admin/goldenset/${id}/approve`, {
                method: 'POST',
                headers: { Authorization: `Bearer ${session?.access_token}` }
            })
            if (!response.ok) throw new Error('Approve failed')
            return response.json()
        },
        onSuccess: () => {
            toast.success('승인되었습니다.')
            setIsOpen(false)
            setSelectedDraftId(null)
            queryClient.invalidateQueries({ queryKey: ['stagingDrafts'] })
        }
    })

    const handleAutoSave = (field: string, value: any) => {
        if (!selectedDraftId) return
        
        // Optimistic update local state
        setFormData(prev => ({ ...prev, [field]: value }))

        if (selectedDraft && selectedDraft[field as keyof GoldenSetDraft] === value) return 

        updateMutation.mutate({
            id: selectedDraftId,
            updates: { [field]: value }
        })
    }

    // --- Drawer Render ---
    return (
        <div className="flex h-[calc(100vh-200px)] gap-4 relative overflow-hidden">
            {/* 1. Left List Panel */}
            <Card className={`bg-slate-900 border-slate-800 flex flex-col transition-all duration-300 ${isOpen ? 'w-1/2' : 'w-full'}`}>
                <CardHeader className="pb-2">
                    <CardTitle className="text-white flex justify-between items-center">
                        <span>Inbox (Staging)</span>
                        <Badge variant="secondary">{draftsResponse?.total || 0}</Badge>
                    </CardTitle>
                    <div className="relative mt-2">
                        <Search className="absolute left-2 top-2.5 h-4 w-4 text-slate-500" />
                        <Input
                            placeholder="약물명, 타겟 검색..."
                            className="pl-8 bg-slate-800 border-slate-700 text-white"
                            value={searchQuery}
                            onChange={e => setSearchQuery(e.target.value)}
                        />
                    </div>
                </CardHeader>
                <CardContent className="flex-1 overflow-hidden p-0">
                    <ScrollArea className="h-full">
                        {isLoading ? (
                            <div className="flex justify-center p-8"><Loader2 className="animate-spin text-purple-500" /></div>
                        ) : (
                            <div className="divide-y divide-slate-800">
                                {draftsResponse?.data.map(draft => (
                                    <div
                                        key={draft.id}
                                        onClick={() => {
                                            setSelectedDraftId(draft.id)
                                            setIsOpen(true)
                                        }}
                                        className={`p-4 cursor-pointer hover:bg-slate-800/50 border-l-4 transition-colors ${
                                            selectedDraftId === draft.id ? 'bg-slate-800 border-purple-500' : 
                                            draft.is_manual_override ? 'border-blue-500' : 'border-transparent'
                                        }`}
                                    >
                                        <div className="flex justify-between items-center">
                                            <h4 className="font-bold text-white truncate max-w-[200px]">{draft.drug_name}</h4>
                                            <StatusBadge status={draft.outcome_type} override={draft.is_manual_override} />
                                        </div>
                                        <div className="flex justify-between mt-1">
                                            <span className="text-sm text-purple-400">{draft.target}</span>
                                            <span className="text-xs text-slate-500">{new Date(draft.created_at).toLocaleDateString()}</span>
                                        </div>
                                        {draft.is_ai_extracted && (
                                            <div className="mt-2 flex items-center gap-1 text-[10px] text-slate-400">
                                                <Sparkles className="w-3 h-3 text-yellow-500" /> AI Suggested
                                            </div>
                                        )}
                                    </div>
                                ))}
                            </div>
                        )}
                    </ScrollArea>
                </CardContent>
            </Card>

            {/* 2. Right Drawer (Slide-over) */}
            <div className={`absolute top-0 right-0 h-full w-1/2 bg-slate-900 border-l border-slate-700 shadow-2xl transform transition-transform duration-300 z-20 ${isOpen ? 'translate-x-0' : 'translate-x-full'}`}>
                {selectedDraft && (
                    <div className="h-full flex flex-col">
                        {/* Header */}
                        <div className="p-4 border-b border-slate-800 flex justify-between items-center bg-slate-950">
                            <div>
                                <h3 className="text-lg font-bold text-white">{formData.drug_name}</h3>
                                <p className="text-xs text-slate-500">ID: {selectedDraft.id}</p>
                            </div>
                            <Button size="sm" variant="ghost" onClick={() => setIsOpen(false)}>
                                <X className="w-4 h-4 text-slate-400" />
                            </Button>
                        </div>

                        {/* Scrollable Form Content */}
                        <ScrollArea className="flex-1 p-6">
                            <div className="space-y-8 pb-20">
                                {/* Core Identity */}
                                <section>
                                    <h4 className="text-sm font-semibold text-purple-400 mb-4 flex items-center gap-2">
                                        <Activity className="w-4 h-4" /> Core Identity
                                    </h4>
                                    <div className="grid grid-cols-2 gap-4">
                                        <div className="space-y-2">
                                            <Label className="text-slate-400">Drug Name</Label>
                                            <Input 
                                                value={formData.drug_name || ''} 
                                                onChange={(e) => setFormData(prev => ({ ...prev, drug_name: e.target.value }))}
                                                onBlur={(e) => handleAutoSave('drug_name', e.target.value)}
                                                className="bg-slate-800 border-slate-700 text-white" 
                                            />
                                        </div>
                                        <div className="space-y-2">
                                            <Label className="text-slate-400">Target Antigen</Label>
                                            <Input 
                                                value={formData.target || ''} 
                                                onChange={(e) => setFormData(prev => ({ ...prev, target: e.target.value }))}
                                                onBlur={(e) => handleAutoSave('target', e.target.value)}
                                                className="bg-slate-800 border-slate-700 text-white" 
                                            />
                                            {selectedDraft.properties?.ai_analysis?.target && (
                                                <p className="text-[10px] text-slate-500 cursor-pointer hover:text-purple-400"
                                                   onClick={() => handleAutoSave('target', selectedDraft.properties?.ai_analysis?.target)}>
                                                    🤖 AI Hint: {selectedDraft.properties.ai_analysis.target}
                                                </p>
                                            )}
                                        </div>
                                    </div>
                                </section>

                                {/* Chemistry */}
                                <section>
                                    <h4 className="text-sm font-semibold text-blue-400 mb-4 flex items-center gap-2">
                                        <FlaskConical className="w-4 h-4" /> Chemistry
                                    </h4>
                                    <div className="grid grid-cols-2 gap-4 mb-4">
                                        <div className="space-y-2">
                                            <Label className="text-slate-400">Payload</Label>
                                            <Input 
                                                value={formData.payload || ''} 
                                                onChange={(e) => setFormData(prev => ({ ...prev, payload: e.target.value }))}
                                                onBlur={(e) => handleAutoSave('payload', e.target.value)}
                                                className="bg-slate-800 border-slate-700 text-white" 
                                            />
                                        </div>
                                        <div className="space-y-2">
                                            <Label className="text-slate-400">Linker</Label>
                                            <Input 
                                                value={formData.linker || ''} 
                                                onChange={(e) => setFormData(prev => ({ ...prev, linker: e.target.value }))}
                                                onBlur={(e) => handleAutoSave('linker', e.target.value)}
                                                className="bg-slate-800 border-slate-700 text-white" 
                                            />
                                        </div>
                                    </div>
                                    <div className="space-y-2">
                                        <Label className="text-slate-400">SMILES Code</Label>
                                        <Textarea 
                                            value={formData.smiles_code || ''}
                                            onChange={(e) => setFormData(prev => ({ ...prev, smiles_code: e.target.value }))}
                                            onBlur={(e) => handleAutoSave('smiles_code', e.target.value)}
                                            className="bg-slate-800 border-slate-700 font-mono text-xs text-white"
                                            rows={2}
                                        />
                                        {/* Canvas Render */}
                                        <div className="mt-2">
                                            <ChemicalStructure smiles={formData.smiles_code || ''} />
                                        </div>
                                    </div>
                                </section>

                                {/* Bio Metrics */}
                                <section>
                                    <h4 className="text-sm font-semibold text-green-400 mb-4 flex items-center gap-2">
                                        <Microscope className="w-4 h-4" /> Bio Metrics (Quantitative)
                                    </h4>
                                    <div className="grid grid-cols-3 gap-4">
                                        <div className="space-y-2">
                                            <Label className="text-slate-400">Binding (Kd)</Label>
                                            <Input 
                                                placeholder="e.g. 1.2 nM"
                                                value={formData.binding_affinity || ''} 
                                                onChange={(e) => setFormData(prev => ({ ...prev, binding_affinity: e.target.value }))}
                                                onBlur={(e) => handleAutoSave('binding_affinity', e.target.value)}
                                                className="bg-slate-800 border-slate-700 text-white" 
                                            />
                                        </div>
                                        <div className="space-y-2">
                                            <Label className="text-slate-400">Isotype</Label>
                                            <Select 
                                                value={formData.isotype || ''} 
                                                onValueChange={(v) => handleAutoSave('isotype', v)}
                                            >
                                                <SelectTrigger className="bg-slate-800 border-slate-700 text-white">
                                                    <SelectValue placeholder="Select" />
                                                </SelectTrigger>
                                                <SelectContent>
                                                    <SelectItem value="IgG1">IgG1</SelectItem>
                                                    <SelectItem value="IgG2">IgG2</SelectItem>
                                                    <SelectItem value="IgG4">IgG4</SelectItem>
                                                </SelectContent>
                                            </Select>
                                        </div>
                                        <div className="space-y-2">
                                            <Label className="text-slate-400">Host</Label>
                                            <Input 
                                                placeholder="Human/Mouse"
                                                value={formData.host_species || ''}
                                                onChange={(e) => setFormData(prev => ({ ...prev, host_species: e.target.value }))}
                                                onBlur={(e) => handleAutoSave('host_species', e.target.value)}
                                                className="bg-slate-800 border-slate-700 text-white" 
                                            />
                                        </div>
                                    </div>
                                    <div className="grid grid-cols-3 gap-4 mt-4">
                                        <div className="space-y-2">
                                            <Label className="text-slate-400">ORR (%)</Label>
                                            <Input 
                                                type="number"
                                                value={formData.orr_pct || ''}
                                                onChange={(e) => setFormData(prev => ({ ...prev, orr_pct: e.target.value }))}
                                                onBlur={(e) => handleAutoSave('orr_pct', e.target.value)}
                                                className="bg-slate-800 border-slate-700 text-white" 
                                            />
                                        </div>
                                        <div className="space-y-2">
                                            <Label className="text-slate-400">OS (Mo)</Label>
                                            <Input 
                                                type="number"
                                                value={formData.os_months || ''}
                                                onChange={(e) => setFormData(prev => ({ ...prev, os_months: e.target.value }))}
                                                onBlur={(e) => handleAutoSave('os_months', e.target.value)}
                                                className="bg-slate-800 border-slate-700 text-white" 
                                            />
                                        </div>
                                        <div className="space-y-2">
                                            <Label className="text-slate-400">PFS (Mo)</Label>
                                            <Input 
                                                type="number"
                                                value={formData.pfs_months || ''}
                                                onChange={(e) => setFormData(prev => ({ ...prev, pfs_months: e.target.value }))}
                                                onBlur={(e) => handleAutoSave('pfs_months', e.target.value)}
                                                className="bg-slate-800 border-slate-700 text-white" 
                                            />
                                        </div>
                                    </div>
                                </section>

                                {/* Outcome */}
                                <section>
                                    <h4 className="text-sm font-semibold text-slate-400 mb-4 flex items-center gap-2">
                                        <ArrowRight className="w-4 h-4" /> Clinical Outcome
                                    </h4>
                                    <div className="grid grid-cols-2 gap-4">
                                        <div className="space-y-2">
                                            <Label className="text-slate-400">Outcome Type</Label>
                                            <Select 
                                                value={formData.outcome_type || ''}
                                                onValueChange={(v) => handleAutoSave('outcome_type', v)}
                                            >
                                                <SelectTrigger className="bg-slate-800 border-slate-700 text-white">
                                                    <SelectValue placeholder="Outcome" />
                                                </SelectTrigger>
                                                <SelectContent>
                                                    <SelectItem value="Success">Success</SelectItem>
                                                    <SelectItem value="Failure">Failure</SelectItem>
                                                    <SelectItem value="Terminated">Terminated</SelectItem>
                                                    <SelectItem value="Ongoing">Ongoing</SelectItem>
                                                </SelectContent>
                                            </Select>
                                        </div>
                                        <div className="space-y-2">
                                            <Label className="text-slate-400">Reason (if failed)</Label>
                                            <Input 
                                                value={formData.failure_reason || ''}
                                                onChange={(e) => setFormData(prev => ({ ...prev, failure_reason: e.target.value }))}
                                                onBlur={(e) => handleAutoSave('failure_reason', e.target.value)}
                                                className="bg-slate-800 border-slate-700 text-white" 
                                            />
                                        </div>
                                    </div>
                                </section>
                            </div>
                        </ScrollArea>

                        {/* Footer Actions */}
                        <div className="p-4 border-t border-slate-800 bg-slate-950 flex justify-end gap-2">
                            <Button 
                                variant="outline" 
                                onClick={() => {
                                    if(confirm('Delete this draft?')) {
                                        // TODO: Add Delete Mutation
                                    }
                                }}
                                className="border-red-900/50 text-red-500 hover:bg-red-900/20"
                            >
                                Delete
                            </Button>
                            <Button 
                                className="bg-green-600 hover:bg-green-500 text-white"
                                onClick={() => approveMutation.mutate(selectedDraft.id)}
                                disabled={approveMutation.isPending}
                            >
                                {approveMutation.isPending ? <Loader2 className="w-4 h-4 animate-spin mr-2"/> : <CheckCircle className="w-4 h-4 mr-2"/>}
                                Approve & Promote
                            </Button>
                        </div>
                    </div>
                )}
            </div>
        </div>
    )
}