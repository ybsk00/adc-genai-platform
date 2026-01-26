import { useState, useEffect } from 'react'
import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Input } from '@/components/ui/input'
import { Badge } from '@/components/ui/badge'
import { Button } from '@/components/ui/button'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Textarea } from '@/components/ui/textarea'
import { Label } from '@/components/ui/label'
import {
    Search, BookOpen, ExternalLink, Loader2, FileText,
    CheckCircle, AlertTriangle, Sparkles, MessageSquare,
    Database, Flame, Wand2, Save
} from 'lucide-react'
import { getSession } from '@/lib/supabase'
import { API_BASE_URL } from '@/lib/api'
import { toast } from 'sonner'
import { AIAssistantPanel } from './AIAssistantPanel'

// --- Types ---
interface KnowledgeItem {
    id: number
    source_type: string
    title: string
    content?: string
    summary?: string
    ai_reasoning?: string
    relevance_score: number
    source_tier?: number
    rag_status: string
    created_at: string
    properties?: any
}

interface KnowledgeEvidenceTabProps {
    onNctClick: (nctId: string) => void
}

export function KnowledgeEvidenceTab({ onNctClick }: KnowledgeEvidenceTabProps) {
    const queryClient = useQueryClient()
    const [searchQuery, setSearchQuery] = useState('')
    const [selectedItemId, setSelectedItemId] = useState<number | null>(null)

    // AI Panel State
    const [isAIPanelOpen, setIsAIPanelOpen] = useState(false)
    const [aiPanelKey, setAiPanelKey] = useState(0)
    const [aiTrigger, setAiTrigger] = useState<string | undefined>(undefined)

    // Edit State
    const [formData, setFormData] = useState<Partial<KnowledgeItem>>({})
    const [isSaving, setIsSaving] = useState(false)
    const [showRawData, setShowRawData] = useState(false)

    // --- Queries ---
    const { data: items, isLoading } = useQuery({
        queryKey: ['knowledgeBase', searchQuery],
        queryFn: async () => {
            const { session } = await getSession()
            const response = await fetch(`${API_BASE_URL}/api/knowledge-base/?limit=50`, {
                headers: { Authorization: `Bearer ${session?.access_token}` }
            })
            if (!response.ok) throw new Error('Failed to fetch knowledge base')
            return response.json() as Promise<KnowledgeItem[]>
        }
    })

    const selectedItem = items?.find(item => item.id === selectedItemId)

    useEffect(() => {
        if (selectedItem) {
            setFormData(selectedItem)
        }
    }, [selectedItem])

    // --- Mutations ---
    const updateMutation = useMutation({
        mutationFn: async ({ id, updates }: { id: number, updates: any }) => {
            const { session } = await getSession()
            const response = await fetch(`${API_BASE_URL}/api/knowledge-base/${id}`, {
                method: 'PATCH',
                headers: {
                    'Content-Type': 'application/json',
                    Authorization: `Bearer ${session?.access_token}`
                },
                body: JSON.stringify(updates)
            })
            if (!response.ok) throw new Error('Update failed')
            return response.json()
        },
        onSuccess: () => {
            toast.success('Saved successfully')
            queryClient.invalidateQueries({ queryKey: ['knowledgeBase'] })
            setIsSaving(false)
        },
        onError: () => {
            toast.error('Save failed')
            setIsSaving(false)
        }
    })

    const handleSave = () => {
        if (!selectedItemId) return
        setIsSaving(true)
        updateMutation.mutate({ id: selectedItemId, updates: formData })
    }

    const toggleAI = () => {
        if (!isAIPanelOpen) {
            setAiPanelKey(prev => prev + 1)
            setIsAIPanelOpen(true)
        } else {
            setIsAIPanelOpen(false)
        }
    }

    const handleMagicAsk = (context: string) => {
        setAiTrigger(context)
        if (!isAIPanelOpen) toggleAI()
    }

    // --- Render Helpers ---
    const renderNctLinks = (text: string) => {
        if (!text) return null
        const nctRegex = /NCT\d{8}/g
        const parts = text.split(nctRegex)
        const matches = text.match(nctRegex) || []

        if (matches.length === 0) return text

        return parts.reduce((acc: React.ReactNode[], part, i) => {
            acc.push(part)
            if (matches[i]) {
                acc.push(
                    <Button
                        key={i}
                        variant="link"
                        className="p-0 h-auto text-blue-400 hover:text-blue-300 font-mono mx-1"
                        onClick={(e) => {
                            e.stopPropagation()
                            onNctClick(matches[i])
                        }}
                    >
                        {matches[i]}
                    </Button>
                )
            }
            return acc
        }, [])
    }

    const filteredItems = items?.filter(item =>
        item.title.toLowerCase().includes(searchQuery.toLowerCase()) ||
        item.summary?.toLowerCase().includes(searchQuery.toLowerCase())
    )

    // --- Main Render ---
    return (
        <div className="flex h-[calc(100vh-250px)] gap-4">
            {/* 1. Left Panel: Master List (35%) */}
            <div className="w-[35%] flex flex-col bg-slate-900 border border-slate-800 rounded-lg overflow-hidden">
                <div className="p-4 border-b border-slate-800 bg-slate-950">
                    <div className="relative">
                        <Search className="absolute left-2 top-2.5 h-4 w-4 text-slate-500" />
                        <Input
                            placeholder="Search title, NCT ID..."
                            className="pl-8 bg-slate-800 border-slate-700 text-white h-9 text-sm"
                            value={searchQuery}
                            onChange={e => setSearchQuery(e.target.value)}
                        />
                    </div>
                    <div className="mt-2 text-xs text-slate-500 flex justify-between">
                        <span>{filteredItems?.length || 0} items found</span>
                        <span>Sort by: Relevance</span>
                    </div>
                </div>

                <ScrollArea className="flex-1">
                    {isLoading ? (
                        <div className="flex justify-center p-8"><Loader2 className="animate-spin text-blue-500" /></div>
                    ) : (
                        <div className="divide-y divide-slate-800">
                            {filteredItems?.map(item => (
                                <div
                                    key={item.id}
                                    onClick={() => setSelectedItemId(item.id)}
                                    className={`p-4 cursor-pointer hover:bg-slate-800/50 transition-colors border-l-4 ${selectedItemId === item.id ? 'bg-slate-800 border-blue-500' : 'border-transparent'
                                        }`}
                                >
                                    <div className="flex justify-between items-start mb-2">
                                        <Badge variant="outline" className="bg-slate-900 text-slate-400 border-slate-700 text-[10px] flex items-center gap-1">
                                            {item.source_type === 'PubMed' ? <BookOpen className="w-3 h-3" /> : <FileText className="w-3 h-3" />}
                                            {item.source_type}
                                        </Badge>
                                        {item.relevance_score >= 0.8 && (
                                            <Badge className="bg-orange-500/20 text-orange-400 border-orange-500/30 text-[10px] flex items-center gap-1">
                                                <Flame className="w-3 h-3" /> Hot Topic
                                            </Badge>
                                        )}
                                    </div>
                                    <h4 className="font-bold text-slate-200 text-sm line-clamp-2 mb-2 leading-snug">
                                        {item.title}
                                    </h4>
                                    <div className="flex justify-between items-center text-xs text-slate-500">
                                        <div className="flex items-center gap-2">
                                            <span className={item.relevance_score >= 0.7 ? "text-green-400" : "text-slate-500"}>
                                                {(item.relevance_score * 100).toFixed(0)}% Match
                                            </span>
                                        </div>
                                        <Badge variant="secondary" className={`text-[10px] h-5 ${item.rag_status === 'indexed' ? 'bg-green-900/20 text-green-400' : 'bg-slate-800 text-slate-500'
                                            }`}>
                                            {item.rag_status}
                                        </Badge>
                                    </div>
                                </div>
                            ))}
                        </div>
                    )}
                </ScrollArea>

                {/* Pagination Placeholder */}
                <div className="p-2 border-t border-slate-800 bg-slate-950 flex justify-center gap-2">
                    <Button variant="ghost" size="sm" disabled>&lt;</Button>
                    <span className="text-xs text-slate-500 self-center">Page 1</span>
                    <Button variant="ghost" size="sm" disabled>&gt;</Button>
                </div>
            </div>

            {/* 2. Right Panel: Detail View (65%) */}
            <div className="flex-1 bg-slate-900 border border-slate-800 rounded-lg overflow-hidden flex flex-col relative">
                {selectedItem ? (
                    <div className="flex h-full">
                        <div className="flex-1 flex flex-col min-w-0">
                            {/* Header */}
                            <div className="p-4 border-b border-slate-800 bg-slate-950 flex justify-between items-start gap-4">
                                <div>
                                    <h3 className="text-lg font-bold text-white leading-tight">
                                        {formData.title}
                                    </h3>
                                    <div className="flex items-center gap-3 mt-2 text-xs text-slate-500">
                                        <span className="flex items-center gap-1">
                                            <BookOpen className="w-3 h-3" /> {formData.source_type}
                                        </span>
                                        <span>•</span>
                                        <span>{new Date(formData.created_at || '').toLocaleDateString()}</span>
                                        <span>•</span>
                                        <span className="font-mono">ID: {formData.id}</span>
                                    </div>
                                </div>
                                <Button size="sm" onClick={handleSave} disabled={isSaving} className="shrink-0 bg-blue-600 hover:bg-blue-500">
                                    {isSaving ? <Loader2 className="w-4 h-4 animate-spin mr-2" /> : <Save className="w-4 h-4 mr-2" />}
                                    Save Changes
                                </Button>
                            </div>

                            <ScrollArea className="flex-1 p-6">
                                <div className="space-y-8 max-w-4xl mx-auto">
                                    {/* Section A: Metadata */}
                                    <section className="space-y-4">
                                        <h4 className="text-sm font-semibold text-slate-400 uppercase tracking-wider border-b border-slate-800 pb-2">
                                            Section A: Metadata
                                        </h4>
                                        <div className="grid grid-cols-3 gap-4">
                                            <div className="space-y-2">
                                                <Label className="text-slate-400 text-xs">Source Type</Label>
                                                <Input value={formData.source_type || ''} readOnly className="bg-slate-950 border-slate-800" />
                                            </div>
                                            <div className="space-y-2">
                                                <Label className="text-slate-400 text-xs">Source Tier</Label>
                                                <Input value={formData.source_tier || '2'} readOnly className="bg-slate-950 border-slate-800" />
                                            </div>
                                            <div className="space-y-2">
                                                <Label className="text-slate-400 text-xs">Relevance Score</Label>
                                                <Input value={formData.relevance_score || 0} readOnly className="bg-slate-950 border-slate-800 font-mono text-blue-400" />
                                            </div>
                                        </div>

                                        {/* Extracted NCT IDs from properties */}
                                        {formData.properties?.nct_id && (
                                            <div className="p-3 bg-blue-900/10 border border-blue-900/30 rounded-lg flex items-center gap-3">
                                                <span className="text-xs font-semibold text-blue-400">Linked Clinical Trial:</span>
                                                <Button
                                                    variant="link"
                                                    className="h-auto p-0 text-blue-300 font-mono text-sm"
                                                    onClick={() => onNctClick(formData.properties.nct_id)}
                                                >
                                                    {formData.properties.nct_id} <ExternalLink className="w-3 h-3 ml-1" />
                                                </Button>
                                            </div>
                                        )}
                                    </section>

                                    {/* Section B: AI Analysis */}
                                    <section className="space-y-4">
                                        <h4 className="text-sm font-semibold text-purple-400 uppercase tracking-wider border-b border-purple-900/30 pb-2 flex items-center gap-2">
                                            <Sparkles className="w-4 h-4" /> Section B: AI Analysis
                                        </h4>

                                        <div className="space-y-2">
                                            <Label className="text-slate-400 text-xs flex justify-between">
                                                Summary
                                                <Button
                                                    variant="ghost" size="sm" className="h-4 p-0 text-purple-400 hover:text-purple-300"
                                                    onClick={() => handleMagicAsk("Summarize this article focusing on ADC efficacy")}
                                                >
                                                    <Wand2 className="w-3 h-3 mr-1" /> Ask AI to Summarize
                                                </Button>
                                            </Label>
                                            <Textarea
                                                value={formData.summary || ''}
                                                onChange={e => setFormData({ ...formData, summary: e.target.value })}
                                                className="bg-slate-950 border-slate-800 min-h-[100px] text-sm leading-relaxed"
                                                placeholder="AI generated summary..."
                                            />
                                        </div>

                                        <div className="space-y-2">
                                            <Label className="text-slate-400 text-xs flex justify-between">
                                                AI Reasoning (Why is this important?)
                                                <Button
                                                    variant="ghost" size="sm" className="h-4 p-0 text-purple-400 hover:text-purple-300"
                                                    onClick={() => handleMagicAsk("Why is this article relevant to our ADC project?")}
                                                >
                                                    <Wand2 className="w-3 h-3 mr-1" /> Ask AI Reasoning
                                                </Button>
                                            </Label>
                                            <Textarea
                                                value={formData.ai_reasoning || ''}
                                                onChange={e => setFormData({ ...formData, ai_reasoning: e.target.value })}
                                                className="bg-slate-950 border-slate-800 min-h-[80px] text-sm"
                                                placeholder="Reasoning..."
                                            />
                                        </div>
                                    </section>

                                    {/* Section C: Full Text */}
                                    <section className="space-y-4">
                                        <h4 className="text-sm font-semibold text-green-400 uppercase tracking-wider border-b border-green-900/30 pb-2 flex items-center gap-2">
                                            <FileText className="w-4 h-4" /> Section C: Full Text & Technical
                                        </h4>

                                        <div className="space-y-2">
                                            <Label className="text-slate-400 text-xs">Content</Label>
                                            <div className="relative">
                                                <Textarea
                                                    value={formData.content || ''}
                                                    onChange={e => setFormData({ ...formData, content: e.target.value })}
                                                    className="bg-slate-950 border-slate-800 min-h-[300px] font-mono text-xs leading-relaxed p-4"
                                                />
                                                <div className="absolute top-2 right-2">
                                                    <Badge variant="outline" className="bg-slate-900 text-slate-500">
                                                        {formData.content?.length || 0} chars
                                                    </Badge>
                                                </div>
                                            </div>
                                        </div>

                                        {/* Raw Data Toggle */}
                                        <div className="pt-4 border-t border-slate-800">
                                            <Button
                                                variant="ghost"
                                                size="sm"
                                                onClick={() => setShowRawData(!showRawData)}
                                                className="text-slate-500 hover:text-white text-xs"
                                            >
                                                <Database className="w-3 h-3 mr-2" />
                                                {showRawData ? "Hide Raw JSON" : "Show Raw JSON"}
                                            </Button>
                                            {showRawData && (
                                                <div className="mt-4 p-4 bg-slate-950 rounded-lg border border-slate-800 overflow-auto max-h-60">
                                                    <pre className="text-[10px] text-slate-400 font-mono">
                                                        {JSON.stringify(formData.properties || {}, null, 2)}
                                                    </pre>
                                                </div>
                                            )}
                                        </div>
                                    </section>
                                </div>
                            </ScrollArea>
                        </div>

                        {/* AI Panel Slide-out */}
                        <div className={`absolute top-0 right-0 h-full w-[400px] bg-slate-950 border-l border-slate-800 shadow-2xl transform transition-transform duration-300 z-20 ${isAIPanelOpen ? 'translate-x-0' : 'translate-x-full'}`}>
                            <AIAssistantPanel
                                key={aiPanelKey}
                                contextData={selectedItem}
                                triggerQuery={aiTrigger}
                                onClearTrigger={() => setAiTrigger(undefined)}
                                onClose={() => setIsAIPanelOpen(false)}
                            />
                        </div>

                        {/* Floating AI Toggle Button */}
                        <Button
                            className={`absolute bottom-6 right-6 h-14 w-14 rounded-full shadow-xl z-10 transition-all duration-300 ${isAIPanelOpen ? 'scale-0 opacity-0' : 'scale-100 opacity-100'}`}
                            onClick={toggleAI}
                            style={{ background: 'linear-gradient(135deg, #3b82f6 0%, #8b5cf6 100%)' }}
                        >
                            <MessageSquare className="w-6 h-6 text-white" />
                        </Button>
                    </div>
                ) : (
                    <div className="flex-1 flex flex-col items-center justify-center text-slate-500">
                        <BookOpen className="w-12 h-12 mb-4 opacity-20" />
                        <p>Select an item from the list to view details.</p>
                    </div>
                )}
            </div>
        </div>
    )
}
