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
    Search, FlaskConical, Microscope, Activity, X, Wand2
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
    const [debouncedSmiles, setDebouncedSmiles] = useState(smiles)
    const [renderMode, setRenderMode] = useState<'canvas' | 'image' | 'error'>('canvas')

    // Debounce logic
    useEffect(() => {
        const timer = setTimeout(() => {
            setDebouncedSmiles(smiles)
            setRenderMode('canvas') // Reset to canvas attempt on change
        }, 300)
        return () => clearTimeout(timer)
    }, [smiles])

    useEffect(() => {
        if (!debouncedSmiles || debouncedSmiles.trim() === '') return
        if (renderMode !== 'canvas') return

        const canvas = canvasRef.current
        if (!canvas) return

        const ctx = canvas.getContext('2d')
        if (!ctx) return

        // Clear canvas
        ctx.clearRect(0, 0, canvas.width, canvas.height)

        try {
            const drawer = new SmilesDrawer.Drawer({
                width: 400,
                height: 250,
                compactDrawing: false,
                terminalCarbons: true,
                explicitHydrogens: false,
                bondThickness: 1.0,
                atomMenu: false
            })

            SmilesDrawer.parse(debouncedSmiles, (tree: any) => {
                if (canvasRef.current && renderMode === 'canvas') {
                    try {
                        drawer.draw(tree, canvasRef.current, 'dark', false)
                    } catch (drawErr) {
                        console.warn('SmilesDrawer Draw Error, switching to fallback:', drawErr)
                        setRenderMode('image')
                    }
                }
            }, (err: any) => {
                console.warn('SMILES Parse Error, switching to fallback:', err)
                setRenderMode('image')
            })
        } catch (e) {
            console.error('SmilesDrawer Init Error:', e)
            setRenderMode('image')
        }
    }, [debouncedSmiles, renderMode])

    if (!smiles) return (
        <div className="h-[250px] flex flex-col items-center justify-center bg-slate-950/50 text-slate-500 border border-dashed border-slate-800 rounded-lg">
            <FlaskConical className="w-8 h-8 mb-2 opacity-20" />
            <span className="text-xs font-medium">No Structure Info (No SMILES)</span>
        </div>
    )

    return (
        <div className="bg-[#1a1b1e] rounded-lg p-4 border border-slate-800 shadow-inner flex justify-center items-center overflow-hidden h-[250px] relative">
            {renderMode === 'canvas' && (
                <canvas ref={canvasRef} width={400} height={250} className="max-w-full max-h-full transition-transform hover:scale-105 duration-300" />
            )}

            {renderMode === 'image' && (
                <img
                    src={`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(debouncedSmiles)}/PNG?record_type=2d&image_size=large`}
                    alt="Structure"
                    className="max-w-full max-h-full object-contain"
                    onError={() => setRenderMode('error')}
                />
            )}

            {renderMode === 'error' && (
                <div className="flex flex-col items-center text-red-500">
                    <FlaskConical className="w-8 h-8 mb-2 opacity-50" />
                    <span className="text-xs">Invalid SMILES Structure</span>
                </div>
            )}
        </div>
    )
}

// 1.5 AI Chat Sidebar
function AISidebar({
    isOpen,
    onClose,
    recordId,
    rawContext,
    onApply
}: {
    isOpen: boolean,
    onClose: () => void,
    recordId: string,
    rawContext?: any,
    onApply: (field: string, value: any) => void
}) {
    const [messages, setMessages] = useState<{ role: 'user' | 'ai', content: string }[]>([])
    const [input, setInput] = useState('')
    const [isTyping, setIsTyping] = useState(false)
    const scrollRef = useRef<HTMLDivElement>(null)

    useEffect(() => {
        if (isOpen && messages.length === 0) {
            setMessages([{
                role: 'ai',
                content: "Hello! Feel free to ask anything about this data. I can analyze it based on the raw data. \nExample: 'What is the Kd value of this drug?', 'Summarize the clinical results'"
            }])
        }
    }, [isOpen])

    useEffect(() => {
        if (scrollRef.current) {
            scrollRef.current.scrollTop = scrollRef.current.scrollHeight
        }
    }, [messages])

    const handleSend = async () => {
        if (!input.trim() || isTyping) return

        const userMsg = input.trim()
        setMessages(prev => [...prev, { role: 'user', content: userMsg }])
        setInput('')
        setIsTyping(true)

        try {
            const { session } = await getSession()
            const response = await fetch(`${API_BASE_URL}/api/admin/ai/chat`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    Authorization: `Bearer ${session?.access_token}`
                },
                body: JSON.stringify({
                    record_id: recordId,
                    message: userMsg,
                    context: rawContext
                })
            })

            if (!response.ok) throw new Error('Failed to chat')
            const data = await response.json()

            setMessages(prev => [...prev, { role: 'ai', content: data.answer }])
        } catch (e) {
            toast.error('AI Response Failed')
        } finally {
            setIsTyping(false)
        }
    }

    return (
        <div className={`fixed top-0 right-0 h-full w-[400px] bg-slate-900 border-l border-slate-700 shadow-2xl z-[100] transform transition-transform duration-300 flex flex-col ${isOpen ? 'translate-x-0' : 'translate-x-full'}`}>
            <div className="p-4 border-b border-slate-800 flex justify-between items-center bg-slate-950">
                <div className="flex items-center gap-2">
                    <Sparkles className="w-5 h-5 text-purple-400" />
                    <span className="font-bold text-white text-sm">AI Assistant</span>
                </div>
                <Button variant="ghost" size="sm" onClick={onClose}>
                    <X className="w-4 h-4 text-slate-400" />
                </Button>
            </div>

            <ScrollArea className="flex-1 p-4" ref={scrollRef}>
                <div className="space-y-4">
                    {messages.map((m, i) => (
                        <div key={i} className={`flex ${m.role === 'user' ? 'justify-end' : 'justify-start'}`}>
                            <div className={`max-w-[85%] p-3 rounded-2xl text-xs leading-relaxed ${m.role === 'user'
                                ? 'bg-purple-600 text-white rounded-tr-none'
                                : 'bg-slate-800 text-slate-200 rounded-tl-none border border-slate-700'
                                }`}>
                                {m.content}
                            </div>
                        </div>
                    ))}
                    {isTyping && (
                        <div className="flex justify-start">
                            <div className="bg-slate-800 p-3 rounded-2xl rounded-tl-none border border-slate-700">
                                <Loader2 className="w-3 h-3 animate-spin text-purple-400" />
                            </div>
                        </div>
                    )}
                </div>
            </ScrollArea>

            <div className="p-4 border-t border-slate-800 bg-slate-950">
                <div className="flex gap-2">
                    <Input
                        placeholder="Type a message..."
                        value={input}
                        onChange={e => setInput(e.target.value)}
                        onKeyDown={e => e.key === 'Enter' && handleSend()}
                        className="bg-slate-800 border-slate-700 text-white text-xs"
                    />
                    <Button size="sm" className="bg-purple-600 hover:bg-purple-500" onClick={handleSend} disabled={isTyping}>
                        Send
                    </Button>
                </div>
            </div>
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
    const [isChatOpen, setIsChatOpen] = useState(false) // AI Chat Sidebar State

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
            toast.success('Saved (Embedding updated)')
            queryClient.invalidateQueries({ queryKey: ['stagingDrafts'] })
        },
        onError: () => toast.error('Save failed')
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
            toast.success('Approved.')
            setIsOpen(false)
            setSelectedDraftId(null)
            queryClient.invalidateQueries({ queryKey: ['stagingDrafts'] })
        }
    })

    const refineMutation = useMutation({
        mutationFn: async (id: string) => {
            const { session } = await getSession()
            const response = await fetch(`${API_BASE_URL}/api/admin/goldenset/${id}/refine`, {
                method: 'POST',
                headers: { Authorization: `Bearer ${session?.access_token}` }
            })
            if (!response.ok) throw new Error('Refinement failed')
            return response.json()
        },
        onSuccess: (data) => {
            toast.success('AI Analysis completed.')
            if (data.analysis) {
                setFormData(prev => ({ ...prev, ...data.analysis }))
            }
            queryClient.invalidateQueries({ queryKey: ['stagingDrafts'] })
        },
        onError: () => toast.error('AI Analysis failed')
    })

    // AI SMILES Fix Mutation
    const fixSmilesMutation = useMutation({
        mutationFn: async () => {
            if (!selectedDraftId || !formData.drug_name) throw new Error('No drug selected')

            const { session } = await getSession()
            const prompt = `
                The user is trying to visualize the chemical structure of '${formData.drug_name}' (Target: ${formData.target || 'Unknown'}).
                The current SMILES code is invalid or missing.
                
                Please provide the CORRECT, CANONICAL SMILES code for '${formData.drug_name}'.
                - If it is an ADC, provide the SMILES for the Payload-Linker complex if possible, or just the Payload.
                - Do NOT include the antibody part as it is too large.
                - Output ONLY the SMILES string. Do not add any markdown, explanations, or quotes.
            `

            const response = await fetch(`${API_BASE_URL}/api/admin/ai/chat`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    Authorization: `Bearer ${session?.access_token}`
                },
                body: JSON.stringify({
                    record_id: selectedDraftId,
                    message: prompt,
                    context: selectedDraft?.raw_data
                })
            })

            if (!response.ok) throw new Error('AI request failed')
            const data = await response.json()
            return data.answer // Assuming 'answer' contains the AI response
        },
        onSuccess: (smiles) => {
            const cleanSmiles = smiles.trim().replace(/`/g, '').replace(/smiles/gi, '').trim()
            setFormData(prev => ({ ...prev, smiles_code: cleanSmiles }))
            handleAutoSave('smiles_code', cleanSmiles)
            toast.success('AI recovered SMILES code.')
        },
        onError: () => toast.error('SMILES recovery failed')
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
                            placeholder="Search drug name, target..."
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
                                        className={`p-4 cursor-pointer hover:bg-slate-800/50 border-l-4 transition-colors ${selectedDraftId === draft.id ? 'bg-slate-800 border-purple-500' :
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
            <div className={`absolute top-0 right-0 h-full w-full md:w-[60vw] lg:w-[45vw] bg-slate-900 border-l border-slate-700 shadow-2xl transform transition-transform duration-300 z-20 ${isOpen ? 'translate-x-0' : 'translate-x-full'}`}>
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
                                {/* AI Analysis Button (Centered) */}
                                <div className="flex justify-center mb-6">
                                    <Button
                                        size="sm"
                                        variant="outline"
                                        className="h-9 px-6 text-xs gap-2 border-purple-500/50 text-purple-400 hover:bg-purple-500/10 hover:text-purple-300 transition-all shadow-[0_0_15px_rgba(168,85,247,0.2)]"
                                        onClick={() => setIsChatOpen(true)}
                                    >
                                        <Sparkles className="w-4 h-4" />
                                        Analyze with AI
                                    </Button>
                                </div>

                                {/* Core Identity */}
                                <section>
                                    <div className="flex items-center mb-4 gap-2">
                                        <h4 className="text-sm font-semibold text-purple-400 flex items-center gap-2">
                                            <Activity className="w-4 h-4" /> Core Identity
                                        </h4>
                                    </div>
                                    <div className="grid grid-cols-2 gap-4">
                                        <div className="space-y-2">
                                            <Label className="text-slate-400 flex items-center gap-1">
                                                Drug Name
                                                <Sparkles className="w-3 h-3 text-purple-500 cursor-pointer hover:text-purple-400" onClick={() => setIsChatOpen(true)} />
                                            </Label>
                                            <Input
                                                value={formData.drug_name || ''}
                                                onChange={(e) => setFormData(prev => ({ ...prev, drug_name: e.target.value }))}
                                                onBlur={(e) => handleAutoSave('drug_name', e.target.value)}
                                                className="bg-slate-800 border-slate-700 text-white"
                                            />
                                        </div>
                                        <div className="space-y-2">
                                            <Label className="text-slate-400 flex items-center gap-1">
                                                Target Antigen
                                                <Sparkles className="w-3 h-3 text-purple-500 cursor-pointer hover:text-purple-400" onClick={() => setIsChatOpen(true)} />
                                            </Label>
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
                                            <Label className="text-slate-400 flex items-center gap-1">
                                                Payload
                                                <Sparkles className="w-3 h-3 text-purple-500 cursor-pointer hover:text-purple-400" onClick={() => setIsChatOpen(true)} />
                                            </Label>
                                            <Input
                                                value={formData.payload || ''}
                                                onChange={(e) => setFormData(prev => ({ ...prev, payload: e.target.value }))}
                                                onBlur={(e) => handleAutoSave('payload', e.target.value)}
                                                className="bg-slate-800 border-slate-700 text-white"
                                            />
                                        </div>
                                        <div className="space-y-2">
                                            <Label className="text-slate-400 flex items-center gap-1">
                                                Linker
                                                <Sparkles className="w-3 h-3 text-purple-500 cursor-pointer hover:text-purple-400" onClick={() => setIsChatOpen(true)} />
                                            </Label>
                                            <Input
                                                value={formData.linker || ''}
                                                onChange={(e) => setFormData(prev => ({ ...prev, linker: e.target.value }))}
                                                onBlur={(e) => handleAutoSave('linker', e.target.value)}
                                                className="bg-slate-800 border-slate-700 text-white"
                                            />
                                        </div>
                                    </div>
                                    <div className="space-y-2">
                                        <Label className="text-slate-400 flex items-center gap-1">
                                            SMILES Code
                                            <Button
                                                size="sm"
                                                variant="ghost"
                                                className="h-5 px-2 text-[10px] bg-purple-500/10 text-purple-400 hover:bg-purple-500/20 hover:text-purple-300 border border-purple-500/30 ml-2"
                                                onClick={() => fixSmilesMutation.mutate()}
                                                disabled={fixSmilesMutation.isPending}
                                            >
                                                {fixSmilesMutation.isPending ? <Loader2 className="w-3 h-3 animate-spin mr-1" /> : <Wand2 className="w-3 h-3 mr-1" />}
                                                AI Structure Fix
                                            </Button>
                                        </Label>
                                        <Textarea
                                            value={formData.smiles_code || ''}
                                            onChange={(e) => setFormData(prev => ({ ...prev, smiles_code: e.target.value }))}
                                            onBlur={(e) => handleAutoSave('smiles_code', e.target.value)}
                                            className="bg-slate-800 border-slate-700 font-mono text-xs text-white"
                                            rows={2}
                                        />
                                        {/* Canvas Render */}
                                        <div className="mt-4">
                                            <Label className="text-[10px] text-slate-500 mb-2 block uppercase tracking-wider">Molecular Structure Visualization</Label>
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
                                            <Label className="text-slate-400 flex items-center gap-1">
                                                Binding (Kd)
                                                <Sparkles className="w-3 h-3 text-purple-500 cursor-pointer hover:text-purple-400" onClick={() => setIsChatOpen(true)} />
                                            </Label>
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
                                            <Label className="text-slate-400 flex items-center gap-1">
                                                ORR (%)
                                                <Sparkles className="w-3 h-3 text-purple-500 cursor-pointer hover:text-purple-400" onClick={() => setIsChatOpen(true)} />
                                            </Label>
                                            <Input
                                                value={formData.orr_pct || ''}
                                                onChange={(e) => setFormData(prev => ({ ...prev, orr_pct: e.target.value }))}
                                                onBlur={(e) => handleAutoSave('orr_pct', e.target.value)}
                                                className="bg-slate-800 border-slate-700 text-white"
                                            />
                                        </div>
                                        <div className="space-y-2">
                                            <Label className="text-slate-400 flex items-center gap-1">
                                                OS (Mo)
                                                <Sparkles className="w-3 h-3 text-purple-500 cursor-pointer hover:text-purple-400" onClick={() => setIsChatOpen(true)} />
                                            </Label>
                                            <Input
                                                value={formData.os_months || ''}
                                                onChange={(e) => setFormData(prev => ({ ...prev, os_months: e.target.value }))}
                                                onBlur={(e) => handleAutoSave('os_months', e.target.value)}
                                                className="bg-slate-800 border-slate-700 text-white"
                                            />
                                        </div>
                                        <div className="space-y-2">
                                            <Label className="text-slate-400 flex items-center gap-1">
                                                PFS (Mo)
                                                <Sparkles className="w-3 h-3 text-purple-500 cursor-pointer hover:text-purple-400" onClick={() => setIsChatOpen(true)} />
                                            </Label>
                                            <Input
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
                                    if (confirm('Delete this draft?')) {
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
                                {approveMutation.isPending ? <Loader2 className="w-4 h-4 animate-spin mr-2" /> : <CheckCircle className="w-4 h-4 mr-2" />}
                                Approve & Promote
                            </Button>
                        </div>
                    </div>
                )}
            </div>

            {/* AI Assistant Sidebar */}
            {selectedDraft && (
                <AISidebar
                    isOpen={isChatOpen}
                    onClose={() => setIsChatOpen(false)}
                    recordId={selectedDraft.id}
                    rawContext={selectedDraft.raw_data}
                    onApply={(field, value) => handleAutoSave(field, value)}
                />
            )}
        </div>
    )
}
