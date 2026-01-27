import { useState, useRef, useEffect } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Textarea } from '@/components/ui/textarea'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Badge } from '@/components/ui/badge'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import {
    X, Send, Bot, Sparkles, Loader2, Edit3, Save, XCircle,
    FileSearch, FlaskConical, BookOpen, Zap, CheckCircle2, AlertTriangle,
    ExternalLink
} from 'lucide-react'
import { cn } from '@/lib/utils'
import { toast } from 'sonner'

interface InventoryDetailPanelProps {
    item: any
    type: 'antibodies' | 'reagents'
    onClose: () => void
    onUpdate?: (updatedItem: any) => void
}

// Quick Actions 프리셋
const QUICK_ACTIONS = {
    validate_smiles: {
        label: "SMILES 검증",
        icon: CheckCircle2,
        prompt: "이 SMILES 구조식이 화학적으로 유효한지 검증해줘. 문제가 있다면 수정된 SMILES를 제안해줘. 관련 PubChem/ChEMBL 레퍼런스 PMID가 있다면 함께 제시해줘.",
        color: "text-green-400"
    },
    find_literature: {
        label: "문헌 검색",
        icon: BookOpen,
        prompt: "이 화합물/항체에 대한 최신 연구 문헌을 찾아줘. 반드시 PMID 또는 DOI 링크를 포함해서 3개 이상의 학술 논문을 제시해줘.",
        color: "text-blue-400"
    },
    suggest_linker: {
        label: "링커 추천",
        icon: FlaskConical,
        prompt: "이 페이로드/항체에 적합한 ADC 링커를 추천해줘. Cleavable vs Non-cleavable 링커의 장단점과 함께 구체적인 링커 구조(SMILES)와 관련 임상 사례(PMID)를 제시해줘.",
        color: "text-purple-400"
    },
    autofill_smiles: {
        label: "SMILES 자동완성",
        icon: Zap,
        prompt: "이 화합물의 CAS 번호 또는 이름을 기반으로 PubChem에서 정확한 SMILES를 찾아줘. MW(분자량)와 함께 검증 결과를 알려줘.",
        color: "text-yellow-400"
    }
}

export function InventoryDetailPanel({ item, type, onClose, onUpdate }: InventoryDetailPanelProps) {
    const [isEditing, setIsEditing] = useState(false)
    const [editData, setEditData] = useState<any>({})
    const [saving, setSaving] = useState(false)
    const [activeTab, setActiveTab] = useState('details')
    const [suggestedSmiles, setSuggestedSmiles] = useState<string | null>(null)

    useEffect(() => {
        if (item) {
            setEditData({ ...item })
            setIsEditing(false)
            setSuggestedSmiles(null)
        }
    }, [item?.id])

    if (!item) {
        return (
            <Card className="h-full bg-slate-900 border-slate-800 flex items-center justify-center text-slate-500">
                <div className="text-center">
                    <Sparkles className="w-12 h-12 mx-auto mb-4 opacity-20" />
                    <p>Select an item to view details</p>
                </div>
            </Card>
        )
    }

    const handleSave = async () => {
        setSaving(true)
        try {
            const API_BASE = import.meta.env.VITE_API_BASE_URL || ''
            const table = type === 'antibodies' ? 'antibody_library' : 'commercial_reagents'
            const res = await fetch(`${API_BASE}/api/admin/inventory/${table}/${item.id}`, {
                method: 'PATCH',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ updates: editData })
            })

            if (res.ok) {
                toast.success('Successfully saved!')
                setIsEditing(false)
                onUpdate?.(editData)
            } else {
                toast.error('Failed to save')
            }
        } catch (error) {
            console.error("Save error:", error)
            toast.error('Network error')
        } finally {
            setSaving(false)
        }
    }

    const handleAcceptSuggestedSmiles = () => {
        if (suggestedSmiles) {
            setEditData({ ...editData, smiles_code: suggestedSmiles })
            setIsEditing(true)
            toast.info('SMILES applied. Click Save to confirm.')
        }
    }

    return (
        <Card className="h-full bg-slate-900 border-slate-800 flex flex-col overflow-hidden">
            {/* Header */}
            <CardHeader className="flex flex-row items-center justify-between py-3 border-b border-slate-800 bg-slate-900/50">
                <div className="flex-1 min-w-0">
                    <CardTitle className="text-lg font-semibold truncate pr-4">
                        {item.product_name || item.name || item.cat_no || 'Unknown'}
                    </CardTitle>
                    <div className="flex gap-2 mt-1">
                        {item.ai_refined && (
                            <Badge className="bg-green-500/20 text-green-300 text-xs">
                                <Sparkles className="w-3 h-3 mr-1" /> AI Refined
                            </Badge>
                        )}
                        {item.is_manual_override && (
                            <Badge className="bg-purple-500/20 text-purple-300 text-xs">
                                Manual Override
                            </Badge>
                        )}
                    </div>
                </div>
                <div className="flex gap-1">
                    {isEditing ? (
                        <>
                            <Button
                                variant="ghost"
                                size="sm"
                                onClick={() => { setIsEditing(false); setEditData({ ...item }); }}
                                className="text-slate-400 hover:text-white"
                            >
                                <XCircle className="w-4 h-4" />
                            </Button>
                            <Button
                                size="sm"
                                onClick={handleSave}
                                disabled={saving}
                                className="bg-green-600 hover:bg-green-700"
                            >
                                {saving ? <Loader2 className="w-4 h-4 animate-spin" /> : <Save className="w-4 h-4" />}
                            </Button>
                        </>
                    ) : (
                        <Button
                            variant="ghost"
                            size="sm"
                            onClick={() => setIsEditing(true)}
                            className="text-slate-400 hover:text-white"
                        >
                            <Edit3 className="w-4 h-4" />
                        </Button>
                    )}
                    <Button variant="ghost" size="icon" onClick={onClose} className="h-8 w-8 text-slate-400 hover:text-white">
                        <X className="w-4 h-4" />
                    </Button>
                </div>
            </CardHeader>

            {/* Tabs */}
            <Tabs value={activeTab} onValueChange={setActiveTab} className="flex-1 flex flex-col min-h-0">
                <TabsList className="grid w-full grid-cols-3 bg-slate-800 rounded-none border-b border-slate-700">
                    <TabsTrigger value="details" className="data-[state=active]:bg-slate-900">Details</TabsTrigger>
                    <TabsTrigger value="structure" className="data-[state=active]:bg-slate-900">Structure</TabsTrigger>
                    <TabsTrigger value="ai" className="data-[state=active]:bg-slate-900">Ask AI</TabsTrigger>
                </TabsList>

                {/* Details Tab */}
                <TabsContent value="details" className="flex-1 overflow-y-auto p-4 space-y-4 m-0">
                    {type === 'antibodies' ? (
                        <AntibodyDetails item={item} editData={editData} setEditData={setEditData} isEditing={isEditing} />
                    ) : (
                        <ReagentDetails item={item} editData={editData} setEditData={setEditData} isEditing={isEditing} />
                    )}
                </TabsContent>

                {/* Structure Comparison Tab */}
                <TabsContent value="structure" className="flex-1 overflow-y-auto p-4 m-0">
                    <StructureComparisonView
                        currentSmiles={item.smiles_code}
                        suggestedSmiles={suggestedSmiles}
                        onAccept={handleAcceptSuggestedSmiles}
                    />
                </TabsContent>

                {/* AI Assistant Tab */}
                <TabsContent value="ai" className="flex-1 flex flex-col min-h-0 m-0">
                    <AIAssistantPanel
                        item={item}
                        type={type}
                        onSmilesGenerated={(smiles) => {
                            setSuggestedSmiles(smiles)
                            setActiveTab('structure')
                        }}
                    />
                </TabsContent>
            </Tabs>
        </Card>
    )
}

// Antibody Details Component
function AntibodyDetails({ item, editData, setEditData, isEditing }: any) {
    return (
        <div className="space-y-3 text-sm">
            <DetailRow label="Cat No" value={isEditing ? editData.cat_no : item.cat_no} field="cat_no" isEditing={isEditing} editData={editData} setEditData={setEditData} />
            <DetailRow label="Product Name" value={isEditing ? editData.product_name : item.product_name} field="product_name" isEditing={isEditing} editData={editData} setEditData={setEditData} />
            <DetailRow label="Target" value={isEditing ? editData.related_disease : item.related_disease} field="related_disease" isEditing={isEditing} editData={editData} setEditData={setEditData} />
            <DetailRow label="Host Species" value={isEditing ? editData.host_species : item.host_species} field="host_species" isEditing={isEditing} editData={editData} setEditData={setEditData} />
            <DetailRow label="Isotype" value={isEditing ? editData.isotype : item.isotype} field="isotype" isEditing={isEditing} editData={editData} setEditData={setEditData} />
            <DetailRow label="Antibody Format" value={isEditing ? editData.antibody_format : item.antibody_format} field="antibody_format" isEditing={isEditing} editData={editData} setEditData={setEditData} />
            <DetailRow label="Source" value={item.source_name} />
            {item.source_url && (
                <div className="pt-2">
                    <a href={item.source_url} target="_blank" rel="noopener noreferrer" className="text-blue-400 hover:underline flex items-center gap-1 text-xs">
                        <ExternalLink className="w-3 h-3" /> View Source
                    </a>
                </div>
            )}
        </div>
    )
}

// Reagent Details Component
function ReagentDetails({ item, editData, setEditData, isEditing }: any) {
    return (
        <div className="space-y-3 text-sm">
            <DetailRow label="Cat No" value={isEditing ? editData.ambeed_cat_no : item.ambeed_cat_no} field="ambeed_cat_no" isEditing={isEditing} editData={editData} setEditData={setEditData} />
            <DetailRow label="Product Name" value={isEditing ? editData.product_name : item.product_name} field="product_name" isEditing={isEditing} editData={editData} setEditData={setEditData} />
            <DetailRow label="Category" value={isEditing ? editData.category : item.category} field="category" isEditing={isEditing} editData={editData} setEditData={setEditData} />
            <DetailRow label="CAS Number" value={isEditing ? editData.cas_number : item.cas_number} field="cas_number" isEditing={isEditing} editData={editData} setEditData={setEditData} />
            <DetailRow label="Formula" value={isEditing ? editData.formula : item.formula} field="formula" isEditing={isEditing} editData={editData} setEditData={setEditData} />
            <DetailRow label="Molecular Weight" value={isEditing ? editData.molecular_weight : item.molecular_weight} field="molecular_weight" isEditing={isEditing} editData={editData} setEditData={setEditData} />
            <DetailRow label="Target" value={isEditing ? editData.target : item.target} field="target" isEditing={isEditing} editData={editData} setEditData={setEditData} />

            {/* SMILES - Special handling with validation indicator */}
            <div className="grid grid-cols-3 gap-2 items-start">
                <span className="text-slate-500">SMILES</span>
                <div className="col-span-2">
                    {isEditing ? (
                        <Textarea
                            value={editData.smiles_code || ''}
                            onChange={(e) => setEditData({ ...editData, smiles_code: e.target.value })}
                            className="bg-slate-800 border-slate-700 text-slate-200 font-mono text-xs min-h-[60px]"
                        />
                    ) : (
                        <div className="flex items-start gap-2">
                            {item.smiles_code ? (
                                <>
                                    <CheckCircle2 className="w-4 h-4 text-green-500 flex-shrink-0 mt-0.5" />
                                    <span className="font-mono text-xs break-all text-slate-200">{item.smiles_code}</span>
                                </>
                            ) : (
                                <span className="text-red-400 flex items-center gap-1">
                                    <AlertTriangle className="w-4 h-4" /> Missing
                                </span>
                            )}
                        </div>
                    )}
                </div>
            </div>

            <DetailRow label="Source" value={item.source_name} />
            <DetailRow label="MDL Number" value={item.mdl_number} />

            {item.product_url && (
                <div className="pt-2">
                    <a href={item.product_url} target="_blank" rel="noopener noreferrer" className="text-blue-400 hover:underline flex items-center gap-1 text-xs">
                        <ExternalLink className="w-3 h-3" /> View Source
                    </a>
                </div>
            )}
        </div>
    )
}

// Detail Row Component (with edit mode support)
function DetailRow({ label, value, field, isEditing, editData, setEditData }: any) {
    if (!value && !isEditing) return null

    return (
        <div className="grid grid-cols-3 gap-2 items-center">
            <span className="text-slate-500">{label}</span>
            {isEditing && field ? (
                <Input
                    value={editData[field] || ''}
                    onChange={(e) => setEditData({ ...editData, [field]: e.target.value })}
                    className="col-span-2 bg-slate-800 border-slate-700 text-slate-200 h-8 text-sm"
                />
            ) : (
                <span className="col-span-2 text-slate-200">{value || '-'}</span>
            )}
        </div>
    )
}

// Structure Comparison View (Current vs AI Suggested)
function StructureComparisonView({ currentSmiles, suggestedSmiles, onAccept }: {
    currentSmiles?: string,
    suggestedSmiles?: string | null,
    onAccept: () => void
}) {
    return (
        <div className="space-y-4">
            <h3 className="text-sm font-medium text-slate-400 uppercase tracking-wider">Structure Comparison</h3>

            <div className="grid grid-cols-2 gap-4">
                {/* Current Structure */}
                <div className="space-y-2">
                    <div className="text-xs text-slate-500 flex items-center gap-1">
                        Current Structure
                        {currentSmiles && <CheckCircle2 className="w-3 h-3 text-green-500" />}
                    </div>
                    <div className="bg-slate-950 border border-slate-800 rounded-lg p-3 min-h-[120px] flex items-center justify-center">
                        {currentSmiles ? (
                            <div className="text-center">
                                <div className="font-mono text-xs text-slate-300 break-all">{currentSmiles}</div>
                                {/* 실제 구조 이미지 렌더링은 RDKit.js나 외부 서비스 필요 */}
                                <div className="mt-2 text-xs text-slate-500">
                                    [Structure visualization requires RDKit.js]
                                </div>
                            </div>
                        ) : (
                            <div className="text-center text-slate-500">
                                <AlertTriangle className="w-8 h-8 mx-auto mb-2 opacity-50" />
                                <p className="text-xs">No structure available</p>
                            </div>
                        )}
                    </div>
                </div>

                {/* AI Suggested Structure */}
                <div className="space-y-2">
                    <div className="text-xs text-slate-500 flex items-center gap-1">
                        AI Suggested
                        {suggestedSmiles && <Sparkles className="w-3 h-3 text-purple-400" />}
                    </div>
                    <div className="bg-slate-950 border border-purple-500/30 rounded-lg p-3 min-h-[120px] flex items-center justify-center">
                        {suggestedSmiles ? (
                            <div className="text-center">
                                <div className="font-mono text-xs text-purple-300 break-all">{suggestedSmiles}</div>
                                <Button
                                    size="sm"
                                    onClick={onAccept}
                                    className="mt-3 bg-purple-600 hover:bg-purple-700"
                                >
                                    <CheckCircle2 className="w-3 h-3 mr-1" /> Accept & Apply
                                </Button>
                            </div>
                        ) : (
                            <div className="text-center text-slate-500">
                                <Bot className="w-8 h-8 mx-auto mb-2 opacity-50" />
                                <p className="text-xs">Use "SMILES 자동완성" in Ask AI</p>
                            </div>
                        )}
                    </div>
                </div>
            </div>

            {/* Comparison Result */}
            {currentSmiles && suggestedSmiles && currentSmiles !== suggestedSmiles && (
                <div className="bg-yellow-500/10 border border-yellow-500/30 rounded-lg p-3">
                    <div className="flex items-center gap-2 text-yellow-300 text-sm">
                        <AlertTriangle className="w-4 h-4" />
                        Structures are different. Review carefully before accepting.
                    </div>
                </div>
            )}
        </div>
    )
}

// AI Assistant Panel with Quick Actions
function AIAssistantPanel({ item, type, onSmilesGenerated }: {
    item: any,
    type: string,
    onSmilesGenerated?: (smiles: string) => void
}) {
    const [messages, setMessages] = useState<{ role: 'user' | 'assistant', content: string }[]>([])
    const [input, setInput] = useState('')
    const [loading, setLoading] = useState(false)
    const scrollRef = useRef<HTMLDivElement>(null)

    useEffect(() => {
        setMessages([])
    }, [item?.id])

    useEffect(() => {
        if (scrollRef.current) {
            scrollRef.current.scrollTop = scrollRef.current.scrollHeight
        }
    }, [messages])

    const handleSend = async (customPrompt?: string) => {
        const userMsg = customPrompt || input
        if (!userMsg.trim() || loading) return

        setMessages(prev => [...prev, { role: 'user', content: userMsg }])
        if (!customPrompt) setInput('')
        setLoading(true)

        try {
            const API_BASE = import.meta.env.VITE_API_BASE_URL || ''
            const res = await fetch(`${API_BASE}/api/admin/ai/chat`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    record_id: item.id,
                    message: userMsg,
                    context: item,
                    mode: 'rag'
                })
            })

            if (res.ok) {
                const data = await res.json()
                const answer = data.answer || "No response generated."
                setMessages(prev => [...prev, { role: 'assistant', content: answer }])

                // SMILES 자동 추출 (AI 응답에서 SMILES 패턴 감지)
                const smilesMatch = answer.match(/SMILES[:\s]+([A-Za-z0-9@+\-\[\]\(\)\\\/=#$%]+)/i)
                if (smilesMatch && onSmilesGenerated) {
                    onSmilesGenerated(smilesMatch[1])
                }
            } else {
                setMessages(prev => [...prev, { role: 'assistant', content: "Error: Failed to get response." }])
            }
        } catch (error) {
            console.error("AI Chat Error", error)
            setMessages(prev => [...prev, { role: 'assistant', content: "Network error. Please try again." }])
        } finally {
            setLoading(false)
        }
    }

    const handleQuickAction = (action: keyof typeof QUICK_ACTIONS) => {
        const prompt = QUICK_ACTIONS[action].prompt
        handleSend(prompt)
    }

    return (
        <div className="flex-1 flex flex-col min-h-0">
            {/* Quick Actions */}
            <div className="p-3 border-b border-slate-800 bg-slate-900/50">
                <div className="text-xs text-slate-500 mb-2">Quick Actions</div>
                <div className="flex flex-wrap gap-2">
                    {Object.entries(QUICK_ACTIONS).map(([key, action]) => (
                        <Button
                            key={key}
                            variant="outline"
                            size="sm"
                            onClick={() => handleQuickAction(key as keyof typeof QUICK_ACTIONS)}
                            disabled={loading}
                            className={cn(
                                "border-slate-700 bg-slate-800 hover:bg-slate-700 text-xs",
                                action.color
                            )}
                        >
                            <action.icon className="w-3 h-3 mr-1" />
                            {action.label}
                        </Button>
                    ))}
                </div>
            </div>

            {/* Chat Area */}
            <ScrollArea className="flex-1 p-4" ref={scrollRef}>
                <div className="space-y-4">
                    {messages.length === 0 && (
                        <div className="text-center text-slate-500 text-sm py-8">
                            <Bot className="w-8 h-8 mx-auto mb-2 opacity-50" />
                            <p>Ask anything about this {type === 'antibodies' ? 'antibody' : 'reagent'}.</p>
                            <p className="text-xs mt-1">Powered by Gemini 2.5 Flash + RAG</p>
                            <p className="text-xs text-purple-400 mt-2">Quick Actions include PMID citations!</p>
                        </div>
                    )}
                    {messages.map((msg, idx) => (
                        <div key={idx} className={cn("flex gap-3", msg.role === 'user' ? "justify-end" : "justify-start")}>
                            {msg.role === 'assistant' && (
                                <div className="w-6 h-6 rounded-full bg-purple-500/20 flex items-center justify-center flex-shrink-0 mt-1">
                                    <Bot className="w-3 h-3 text-purple-400" />
                                </div>
                            )}
                            <div className={cn(
                                "rounded-lg px-3 py-2 max-w-[85%] text-sm whitespace-pre-wrap",
                                msg.role === 'user'
                                    ? "bg-blue-600 text-white"
                                    : "bg-slate-800 text-slate-200"
                            )}>
                                {msg.content}
                            </div>
                        </div>
                    ))}
                    {loading && (
                        <div className="flex gap-3 justify-start">
                            <div className="w-6 h-6 rounded-full bg-purple-500/20 flex items-center justify-center flex-shrink-0 mt-1">
                                <Bot className="w-3 h-3 text-purple-400" />
                            </div>
                            <div className="bg-slate-800 rounded-lg px-3 py-2">
                                <Loader2 className="w-4 h-4 animate-spin text-slate-400" />
                            </div>
                        </div>
                    )}
                </div>
            </ScrollArea>

            {/* Input Area */}
            <div className="p-3 border-t border-slate-800 flex gap-2">
                <Input
                    value={input}
                    onChange={(e) => setInput(e.target.value)}
                    onKeyDown={(e) => e.key === 'Enter' && !e.shiftKey && handleSend()}
                    placeholder="Ask a question..."
                    className="bg-slate-800 border-slate-700 text-slate-200 focus:ring-purple-500"
                    disabled={loading}
                />
                <Button
                    size="icon"
                    onClick={() => handleSend()}
                    disabled={loading || !input.trim()}
                    className="bg-purple-600 hover:bg-purple-700"
                >
                    <Send className="w-4 h-4" />
                </Button>
            </div>
        </div>
    )
}
