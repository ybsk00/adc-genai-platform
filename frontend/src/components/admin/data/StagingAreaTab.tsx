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
    Search, FlaskConical, Microscope, Activity, X, Wand2,
    FileText, AlertTriangle, Copy, Database
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
    raw_data?: any
    // Bio Metrics
    binding_affinity?: string
    isotype?: string
    host_species?: string
    orr_pct?: string
    os_months?: string
    pfs_months?: string
    dor_months?: string
    patient_count?: string
    adverse_events_grade3_pct?: string
    // ADC Design
    target_1?: string
    target_2?: string
    target_symbol?: string
    antibody_format?: string
    linker_type?: string
    dar?: string
    gene_id?: string
    uniprot_id?: string
    // Chemical / IP
    payload_smiles?: string
    linker_smiles?: string
    full_smiles?: string
    canonical_smiles?: string
    molecular_weight?: string
    patent_id?: string
    patent_expiry?: string
    // Meta
    is_manual_override?: boolean
    properties?: Record<string, any>
    status?: string
    reviewer_note?: string
}

interface DraftsResponse {
    data: GoldenSetDraft[]
    total: number
    limit: number
    offset: number
}

// --- Components ---

// 1. Chemical Structure Renderer (Kept same)
function ChemicalStructure({ smiles }: { smiles: string }) {
    const canvasRef = useRef<HTMLCanvasElement>(null)
    const [debouncedSmiles, setDebouncedSmiles] = useState(smiles)
    const [renderMode, setRenderMode] = useState<'canvas' | 'image' | 'error'>('canvas')
    const [canvasId] = useState(() => `smiles-canvas-${Math.random().toString(36).substr(2, 9)}`)

    useEffect(() => {
        const timer = setTimeout(() => {
            setDebouncedSmiles(smiles)
            setRenderMode('canvas')
        }, 300)
        return () => clearTimeout(timer)
    }, [smiles])

    useEffect(() => {
        if (!debouncedSmiles || debouncedSmiles.trim() === '') return
        if (renderMode !== 'canvas') return

        const canvas = canvasRef.current
        if (!canvas) return

        const renderId = requestAnimationFrame(() => {
            try {
                const ctx = canvas.getContext('2d')
                if (ctx) ctx.clearRect(0, 0, canvas.width, canvas.height)

                const DrawerClass = (SmilesDrawer as any).SmiDrawer || (SmilesDrawer as any).Drawer
                if (!DrawerClass) throw new Error('SmilesDrawer class not found')

                const drawer = new DrawerClass({
                    width: 800, height: 500, compactDrawing: false,
                    terminalCarbons: true, explicitHydrogens: false,
                    bondThickness: 2.0, bondLength: 30, fontSize: 24, theme: 'dark'
                })

                if ((SmilesDrawer as any).SmiDrawer) {
                    drawer.draw(debouncedSmiles, `#${canvasId}`, 'dark')
                } else {
                    SmilesDrawer.parse(debouncedSmiles, (tree: any) => {
                        drawer.draw(tree, canvas, 'dark', false)
                    }, () => setRenderMode('image'))
                }
            } catch (e) {
                setRenderMode('image')
            }
        })
        return () => cancelAnimationFrame(renderId)
    }, [debouncedSmiles, renderMode, canvasId])

    if (!smiles) return (
        <div className="h-[150px] flex flex-col items-center justify-center bg-slate-950/50 text-slate-500 border border-dashed border-slate-800 rounded-lg">
            <FlaskConical className="w-6 h-6 mb-2 opacity-20" />
            <span className="text-xs font-medium">No Structure Info</span>
        </div>
    )

    return (
        <div className="bg-[#1a1b1e] rounded-lg p-2 border border-slate-800 shadow-inner flex justify-center items-center overflow-hidden h-[150px] relative group">
            {renderMode === 'canvas' && (
                <canvas id={canvasId} ref={canvasRef} width={800} height={500} className="max-w-full max-h-full transition-transform group-hover:scale-110 duration-300" />
            )}
            {renderMode === 'image' && (
                <img src={`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(debouncedSmiles)}/PNG?record_type=2d&image_size=large`} alt="Structure" className="max-w-full max-h-full object-contain" onError={() => setRenderMode('error')} />
            )}
            {renderMode === 'error' && (
                <div className="flex flex-col items-center text-red-500">
                    <FlaskConical className="w-6 h-6 mb-2 opacity-50" />
                    <span className="text-xs">Invalid SMILES</span>
                </div>
            )}
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

// 3. Metric Value Display (Visual Cues)
function MetricValue({ value, unit, type }: { value?: string | number, unit?: string, type: 'orr' | 'pfs' | 'os' }) {
    if (!value) return <span className="text-red-500/70 text-xs italic">Missing</span>

    const numVal = parseFloat(String(value))
    let colorClass = "text-slate-200"

    if (type === 'orr') {
        if (numVal >= 50) colorClass = "text-blue-400 font-bold"
        else if (numVal < 30) colorClass = "text-orange-400"
    }

    return <span className={`font-mono ${colorClass}`}>{value}{unit}</span>
}

// ... imports ...

interface StagingAreaTabProps {
    initialSearchQuery?: string
    onSearchClear?: () => void
}

export function StagingAreaTab({ initialSearchQuery, onSearchClear }: StagingAreaTabProps) {
    const queryClient = useQueryClient()
    const [selectedDraftId, setSelectedDraftId] = useState<string | null>(null)
    const [page, setPage] = useState(1)
    const [searchQuery, setSearchQuery] = useState(initialSearchQuery || '')
    const [debouncedSearch, setDebouncedSearch] = useState(initialSearchQuery || '')
    const [showRawData, setShowRawData] = useState(false)

    // Update search when prop changes (from cross-tab link)
    useEffect(() => {
        if (initialSearchQuery) {
            setSearchQuery(initialSearchQuery)
            setDebouncedSearch(initialSearchQuery)
        }
    }, [initialSearchQuery])

    // Form State
    const [formData, setFormData] = useState<Partial<GoldenSetDraft>>({})
    const [isSaving, setIsSaving] = useState(false)

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
                limit: '20', offset: String(offset), search: debouncedSearch
            })
            const response = await fetch(`${API_BASE_URL}/api/admin/goldenset/drafts?${params}`, {
                headers: { Authorization: `Bearer ${session?.access_token}` }
            })
            if (!response.ok) throw new Error('Failed to fetch')
            return response.json() as Promise<DraftsResponse>
        }
    })

    const selectedDraft = draftsResponse?.data.find(d => d.id === selectedDraftId)

    useEffect(() => {
        if (selectedDraft) {
            // Flatten properties into formData for easier editing
            const flattened = { ...selectedDraft, ...selectedDraft.properties }
            setFormData(flattened)
        }
    }, [selectedDraft])

    // --- Mutations ---
    const updateMutation = useMutation({
        mutationFn: async ({ id, updates }: { id: string, updates: any }) => {
            const { session } = await getSession()
            // Use golden_set_library endpoint (via generic inventory patch)
            const response = await fetch(`${API_BASE_URL}/api/admin/inventory/golden_set_library/${id}`, {
                method: 'PATCH',
                headers: { 'Content-Type': 'application/json', Authorization: `Bearer ${session?.access_token}` },
                body: JSON.stringify({ updates })
            })
            if (!response.ok) throw new Error('Update failed')
            return response.json()
        },
        onSuccess: () => {
            toast.success('Saved successfully')
            queryClient.invalidateQueries({ queryKey: ['stagingDrafts'] })
            setIsSaving(false)
        },
        onError: () => {
            toast.error('Save failed')
            setIsSaving(false)
        }
    })

    const handleSave = () => {
        if (!selectedDraftId) return
        setIsSaving(true)
        updateMutation.mutate({ id: selectedDraftId, updates: formData })
    }

    const handlePromote = async () => {
        if (!selectedDraftId) return
        // Just update status to 'approved' via same patch for now, or use specific endpoint if needed
        // Plan says: "Ensure status update to 'approved' moves item"
        // We can use the patch endpoint to set status='approved'
        updateMutation.mutate({
            id: selectedDraftId,
            updates: { ...formData, status: 'approved' }
        })
        toast.success("Promoted to Golden Set!")
        setSelectedDraftId(null) // Close panel
    }

    // --- Render ---
    return (
        <div className="flex h-[calc(100vh-250px)] gap-4">
            {/* 1. Left Panel: Master List */}
            <div className="w-1/3 flex flex-col bg-slate-900 border border-slate-800 rounded-lg overflow-hidden">
                <div className="p-4 border-b border-slate-800 bg-slate-950">
                    <div className="relative">
                        <Search className="absolute left-2 top-2.5 h-4 w-4 text-slate-500" />
                        <Input
                            placeholder="Search drug name, target..."
                            className="pl-8 bg-slate-800 border-slate-700 text-white h-9 text-sm"
                            value={searchQuery}
                            onChange={e => {
                                setSearchQuery(e.target.value)
                                if (e.target.value === '' && onSearchClear) onSearchClear()
                            }}
                        />
                        {searchQuery && (
                            <Button
                                variant="ghost"
                                size="sm"
                                className="absolute right-1 top-1 h-7 w-7 p-0 text-slate-500 hover:text-white"
                                onClick={() => {
                                    setSearchQuery('')
                                    if (onSearchClear) onSearchClear()
                                }}
                            >
                                <X className="h-3 w-3" />
                            </Button>
                        )}
                    </div>
                </div>
                <ScrollArea className="flex-1">
                    {isLoading ? (
                        <div className="flex justify-center p-8"><Loader2 className="animate-spin text-purple-500" /></div>
                    ) : (
                        <div className="divide-y divide-slate-800">
                            {draftsResponse?.data.map(draft => (
                                <div
                                    key={draft.id}
                                    onClick={() => setSelectedDraftId(draft.id)}
                                    className={`p-4 cursor-pointer hover:bg-slate-800/50 transition-colors border-l-4 ${selectedDraftId === draft.id ? 'bg-slate-800 border-purple-500' : 'border-transparent'
                                        }`}
                                >
                                    <div className="flex justify-between items-start mb-1">
                                        <h4 className="font-bold text-slate-200 truncate max-w-[180px] text-sm">{draft.drug_name}</h4>
                                        <StatusBadge status={draft.outcome_type} override={draft.is_manual_override} />
                                    </div>
                                    <div className="flex justify-between items-center text-xs text-slate-500 mb-2">
                                        <span>{draft.target || 'Unknown Target'}</span>
                                        <span>{new Date(draft.created_at).toLocaleDateString()}</span>
                                    </div>
                                    {/* Key Metrics Preview */}
                                    <div className="flex gap-3 text-xs">
                                        <div className="flex items-center gap-1">
                                            <span className="text-slate-600">ORR</span>
                                            <MetricValue value={draft.orr_pct || draft.properties?.orr_pct} unit="%" type="orr" />
                                        </div>
                                        <div className="flex items-center gap-1">
                                            <span className="text-slate-600">PFS</span>
                                            <MetricValue value={draft.pfs_months || draft.properties?.pfs_months} unit="m" type="pfs" />
                                        </div>
                                    </div>
                                </div>
                            ))}
                        </div>
                    )}
                </ScrollArea>
                {/* Pagination */}
                <div className="p-2 border-t border-slate-800 bg-slate-950 flex justify-center gap-2">
                    <Button variant="ghost" size="sm" disabled={page === 1} onClick={() => setPage(p => p - 1)}>&lt;</Button>
                    <span className="text-xs text-slate-500 self-center">Page {page}</span>
                    <Button variant="ghost" size="sm" disabled={!draftsResponse || draftsResponse.data.length < 20} onClick={() => setPage(p => p + 1)}>&gt;</Button>
                </div>
            </div>

            {/* 2. Right Panel: Detail View */}
            <div className="flex-1 bg-slate-900 border border-slate-800 rounded-lg overflow-hidden flex flex-col">
                {selectedDraft ? (
                    <>
                        {/* Header */}
                        <div className="p-4 border-b border-slate-800 bg-slate-950 flex justify-between items-center">
                            <div>
                                <h3 className="text-lg font-bold text-white flex items-center gap-2">
                                    {formData.drug_name}
                                    {formData.is_ai_extracted && <Sparkles className="w-4 h-4 text-yellow-500" />}
                                </h3>
                                <p className="text-xs text-slate-500 font-mono">{selectedDraft.id}</p>
                            </div>
                            <div className="flex gap-2">
                                <Button size="sm" variant="outline" onClick={handleSave} disabled={isSaving} className="h-8 text-xs">
                                    {isSaving ? <Loader2 className="w-3 h-3 animate-spin mr-1" /> : <CheckCircle className="w-3 h-3 mr-1" />}
                                    Save Changes
                                </Button>
                                <Button size="sm" onClick={handlePromote} className="h-8 text-xs bg-purple-600 hover:bg-purple-500">
                                    Promote to Golden Set
                                </Button>
                            </div>
                        </div>

                        <ScrollArea className="flex-1 p-6">
                            <div className="space-y-8 max-w-4xl mx-auto">
                                {/* Section A: Basic Info */}
                                <section className="space-y-4">
                                    <h4 className="text-sm font-semibold text-slate-400 uppercase tracking-wider border-b border-slate-800 pb-2">
                                        Section A: Basic & Clinical Status
                                    </h4>
                                    <div className="grid grid-cols-2 gap-6">
                                        <div className="space-y-2">
                                            <Label className="text-slate-400 text-xs">Drug Name</Label>
                                            <Input
                                                value={formData.drug_name || ''}
                                                onChange={e => setFormData({ ...formData, drug_name: e.target.value })}
                                                className="bg-slate-950 border-slate-800"
                                            />
                                        </div>
                                        <div className="space-y-2">
                                            <Label className="text-slate-400 text-xs">Outcome Type</Label>
                                            <Select
                                                value={formData.outcome_type || 'Unknown'}
                                                onValueChange={v => setFormData({ ...formData, outcome_type: v })}
                                            >
                                                <SelectTrigger className="bg-slate-950 border-slate-800">
                                                    <SelectValue />
                                                </SelectTrigger>
                                                <SelectContent>
                                                    <SelectItem value="Success">Success (Approved)</SelectItem>
                                                    <SelectItem value="Failure">Failure</SelectItem>
                                                    <SelectItem value="Terminated">Terminated</SelectItem>
                                                    <SelectItem value="Unknown">Unknown</SelectItem>
                                                </SelectContent>
                                            </Select>
                                        </div>
                                        <div className="col-span-2 space-y-2">
                                            <Label className="text-slate-400 text-xs">Reviewer Note</Label>
                                            <Textarea
                                                value={formData.reviewer_note || ''}
                                                onChange={e => setFormData({ ...formData, reviewer_note: e.target.value })}
                                                className="bg-slate-950 border-slate-800 h-20 text-xs"
                                                placeholder="Add notes for other reviewers..."
                                            />
                                        </div>
                                    </div>
                                </section>

                                {/* Section B: Metrics */}
                                <section className="space-y-4">
                                    <h4 className="text-sm font-semibold text-blue-400 uppercase tracking-wider border-b border-blue-900/30 pb-2 flex items-center gap-2">
                                        <Activity className="w-4 h-4" /> Section B: Clinical Metrics
                                    </h4>
                                    <div className="grid grid-cols-3 gap-4">
                                        {[
                                            { label: 'ORR (%)', key: 'orr_pct', type: 'orr' },
                                            { label: 'PFS (mo)', key: 'pfs_months', type: 'pfs' },
                                            { label: 'OS (mo)', key: 'os_months', type: 'os' },
                                            { label: 'DoR (mo)', key: 'dor_months', type: 'pfs' },
                                            { label: 'Patient Count', key: 'patient_count', type: 'os' },
                                            { label: 'AE Grade 3+ (%)', key: 'adverse_events_grade3_pct', type: 'orr' },
                                        ].map((field) => (
                                            <div key={field.key} className="space-y-2">
                                                <Label className="text-slate-400 text-xs">{field.label}</Label>
                                                <div className="relative">
                                                    <Input
                                                        value={formData[field.key as keyof GoldenSetDraft] || ''}
                                                        onChange={e => setFormData({ ...formData, [field.key]: e.target.value })}
                                                        className={`bg-slate-950 border-slate-800 ${!formData[field.key as keyof GoldenSetDraft] ? 'border-red-900/50 bg-red-950/10' : ''}`}
                                                        placeholder="-"
                                                    />
                                                    {/* Visual Cue Overlay (Only if not editing? For now just input) */}
                                                </div>
                                            </div>
                                        ))}
                                    </div>
                                </section>

                                {/* Section C: ADC Design */}
                                <section className="space-y-4">
                                    <h4 className="text-sm font-semibold text-purple-400 uppercase tracking-wider border-b border-purple-900/30 pb-2 flex items-center gap-2">
                                        <Microscope className="w-4 h-4" /> Section C: ADC Design
                                    </h4>
                                    <div className="grid grid-cols-3 gap-4">
                                        <div className="space-y-2">
                                            <Label className="text-slate-400 text-xs">Target Antigen</Label>
                                            <Input
                                                value={formData.target || ''}
                                                onChange={e => setFormData({ ...formData, target: e.target.value })}
                                                className="bg-slate-950 border-slate-800"
                                            />
                                        </div>
                                        <div className="space-y-2">
                                            <Label className="text-slate-400 text-xs">Antibody Format</Label>
                                            <Input
                                                value={formData.antibody_format || ''}
                                                onChange={e => setFormData({ ...formData, antibody_format: e.target.value })}
                                                className="bg-slate-950 border-slate-800"
                                                placeholder="e.g. IgG1"
                                            />
                                        </div>
                                        <div className="space-y-2">
                                            <Label className="text-slate-400 text-xs">Linker Type</Label>
                                            <Input
                                                value={formData.linker_type || ''}
                                                onChange={e => setFormData({ ...formData, linker_type: e.target.value })}
                                                className="bg-slate-950 border-slate-800"
                                                placeholder="e.g. Cleavable"
                                            />
                                        </div>
                                        <div className="space-y-2">
                                            <Label className="text-slate-400 text-xs">DAR</Label>
                                            <Input
                                                value={formData.dar || ''}
                                                onChange={e => setFormData({ ...formData, dar: e.target.value })}
                                                className="bg-slate-950 border-slate-800"
                                            />
                                        </div>
                                        <div className="space-y-2">
                                            <Label className="text-slate-400 text-xs">Binding (Kd)</Label>
                                            <Input
                                                value={formData.binding_affinity || ''}
                                                onChange={e => setFormData({ ...formData, binding_affinity: e.target.value })}
                                                className="bg-slate-950 border-slate-800"
                                            />
                                        </div>
                                    </div>
                                </section>

                                {/* Section D: Chemical / IP */}
                                <section className="space-y-4">
                                    <h4 className="text-sm font-semibold text-green-400 uppercase tracking-wider border-b border-green-900/30 pb-2 flex items-center gap-2">
                                        <FlaskConical className="w-4 h-4" /> Section D: Chemical & IP
                                    </h4>
                                    <div className="grid grid-cols-2 gap-6">
                                        <div className="space-y-2">
                                            <Label className="text-slate-400 text-xs flex justify-between">
                                                SMILES Code
                                                <Button variant="ghost" size="sm" className="h-4 p-0 text-slate-500 hover:text-white" onClick={() => {
                                                    navigator.clipboard.writeText(formData.smiles_code || '')
                                                    toast.success("Copied SMILES")
                                                }}>
                                                    <Copy className="w-3 h-3 mr-1" /> Copy
                                                </Button>
                                            </Label>
                                            <Textarea
                                                value={formData.smiles_code || ''}
                                                onChange={e => setFormData({ ...formData, smiles_code: e.target.value })}
                                                className="bg-slate-950 border-slate-800 font-mono text-xs h-20"
                                            />
                                            <ChemicalStructure smiles={formData.smiles_code || ''} />
                                        </div>
                                        <div className="space-y-4">
                                            <div className="space-y-2">
                                                <Label className="text-slate-400 text-xs">Molecular Weight</Label>
                                                <Input
                                                    value={formData.molecular_weight || ''}
                                                    onChange={e => setFormData({ ...formData, molecular_weight: e.target.value })}
                                                    className="bg-slate-950 border-slate-800"
                                                />
                                            </div>
                                            <div className="space-y-2">
                                                <Label className="text-slate-400 text-xs">Patent ID</Label>
                                                <Input
                                                    value={formData.patent_id || ''}
                                                    onChange={e => setFormData({ ...formData, patent_id: e.target.value })}
                                                    className="bg-slate-950 border-slate-800"
                                                />
                                            </div>
                                            <div className="space-y-2">
                                                <Label className="text-slate-400 text-xs">Patent Expiry</Label>
                                                <Input
                                                    value={formData.patent_expiry || ''}
                                                    onChange={e => setFormData({ ...formData, patent_expiry: e.target.value })}
                                                    className="bg-slate-950 border-slate-800"
                                                    type="date"
                                                />
                                            </div>
                                        </div>
                                    </div>
                                </section>

                                {/* Raw Data Toggle */}
                                <div className="pt-8 border-t border-slate-800">
                                    <Button
                                        variant="ghost"
                                        size="sm"
                                        onClick={() => setShowRawData(!showRawData)}
                                        className="text-slate-500 hover:text-white text-xs"
                                    >
                                        <Database className="w-3 h-3 mr-2" />
                                        {showRawData ? "Hide Raw Data" : "Show Raw Data (JSON)"}
                                    </Button>
                                    {showRawData && (
                                        <div className="mt-4 p-4 bg-slate-950 rounded-lg border border-slate-800 overflow-auto max-h-60">
                                            <pre className="text-[10px] text-slate-400 font-mono">
                                                {JSON.stringify(selectedDraft.raw_data || {}, null, 2)}
                                            </pre>
                                        </div>
                                    )}
                                </div>
                            </div>
                        </ScrollArea>
                    </>
                ) : (
                    <div className="flex-1 flex flex-col items-center justify-center text-slate-500">
                        <Microscope className="w-12 h-12 mb-4 opacity-20" />
                        <p>Select a candidate from the list to view details.</p>
                    </div>
                )}
            </div>
        </div>
    )
}
