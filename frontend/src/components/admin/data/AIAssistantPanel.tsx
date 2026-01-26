import { useState, useEffect, useRef } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Switch } from '@/components/ui/switch'
import { Label } from '@/components/ui/label'
import { Bot, User, Sparkles, Database, ExternalLink, Wand2, AlertTriangle, X } from 'lucide-react'
import { API_BASE_URL } from '@/lib/api'
import { getSession } from '@/lib/supabase'
import { toast } from 'sonner'

interface Message {
    role: 'user' | 'assistant'
    content: string
    mode?: 'rag' | 'general'
    sources?: any[]
}

interface AIAssistantPanelProps {
    contextData: any
    onSourceClick?: (sourceType: string, sourceId: string) => void
    triggerQuery?: string // Prop to trigger a query externally (Magic Fill)
    onClearTrigger?: () => void
    onClose?: () => void
}

export function AIAssistantPanel({ contextData, onSourceClick, triggerQuery, onClearTrigger, onClose }: AIAssistantPanelProps) {
    const [messages, setMessages] = useState<Message[]>([])
    const [input, setInput] = useState('')
    const [isLoading, setIsLoading] = useState(false)
    const [mode, setMode] = useState<'rag' | 'general'>('rag')
    const scrollRef = useRef<HTMLDivElement>(null)

    // Handle External Trigger (Magic Fill)
    useEffect(() => {
        if (triggerQuery) {
            handleSend(triggerQuery, 'rag', true) // Force RAG for magic fill
            if (onClearTrigger) onClearTrigger()
        }
    }, [triggerQuery])

    // Auto-scroll
    useEffect(() => {
        if (scrollRef.current) {
            scrollRef.current.scrollTop = scrollRef.current.scrollHeight
        }
    }, [messages])

    const handleSend = async (text: string, activeMode: 'rag' | 'general' = mode, isAutofill = false) => {
        if (!text.trim()) return

        const newMessage: Message = { role: 'user', content: text, mode: activeMode }
        setMessages(prev => [...prev, newMessage])
        setInput('')
        setIsLoading(true)

        try {
            const { session } = await getSession()
            const response = await fetch(`${API_BASE_URL}/api/admin/ai/chat`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    Authorization: `Bearer ${session?.access_token}`
                },
                body: JSON.stringify({
                    record_id: contextData.id || 'unknown',
                    message: text,
                    context: contextData,
                    mode: activeMode,
                    target_field: isAutofill ? 'autofill_target' : undefined // Simplified for now, or pass actual field name if we had it
                })
            })

            if (!response.ok) throw new Error('Failed to get response')

            const data = await response.json()

            setMessages(prev => [...prev, {
                role: 'assistant',
                content: data.answer,
                mode: data.mode,
                sources: data.sources
            }])

        } catch (error) {
            toast.error("AI failed to respond")
            setMessages(prev => [...prev, { role: 'assistant', content: "Sorry, I encountered an error." }])
        } finally {
            setIsLoading(false)
        }
    }

    // Render content with clickable citations
    const renderContent = (msg: Message) => {
        if (msg.role === 'user') return <p className="text-sm text-white">{msg.content}</p>

        // Split by [Source: type:id]
        const parts = msg.content.split(/(\[Source: [^\]]+\])/g)

        return (
            <div className="text-sm text-slate-300 space-y-2">
                <p>
                    {parts.map((part, i) => {
                        const match = part.match(/\[Source: ([^:]+):([^\]]+)\]/)
                        if (match) {
                            const [_, type, id] = match
                            return (
                                <Button
                                    key={i}
                                    variant="outline"
                                    size="sm"
                                    className="h-5 px-1.5 text-[10px] mx-1 bg-blue-900/20 border-blue-800 text-blue-300 hover:bg-blue-800/50"
                                    onClick={() => onSourceClick?.(type, id)}
                                >
                                    <Database className="w-3 h-3 mr-1" />
                                    {type}:{id.substring(0, 6)}...
                                </Button>
                            )
                        }
                        return part
                    })}
                </p>

                {/* Honest AI Fallback */}
                {msg.mode === 'rag' && msg.content.includes("데이터가 없어") && (
                    <div className="mt-2 p-2 bg-yellow-900/20 border border-yellow-800/50 rounded flex items-center gap-2">
                        <AlertTriangle className="w-4 h-4 text-yellow-500" />
                        <span className="text-xs text-yellow-200">No internal data found.</span>
                        <Button
                            variant="link"
                            className="text-yellow-400 h-auto p-0 text-xs underline"
                            onClick={() => handleSend(messages[messages.length - 2].content, 'general')}
                        >
                            Try General Mode?
                        </Button>
                    </div>
                )}
            </div>
        )
    }

    return (
        <Card className="bg-slate-900 border-slate-800 h-full flex flex-col shadow-xl">
            <CardHeader className="py-3 px-4 border-b border-slate-800 flex flex-row items-center justify-between">
                <CardTitle className="text-sm font-medium flex items-center gap-2 text-slate-200">
                    <Bot className="w-4 h-4 text-purple-400" />
                    AI Consultation
                </CardTitle>
                <div className="flex items-center gap-2">
                    <Label htmlFor="mode-switch" className="text-[10px] text-slate-400 font-mono">
                        {mode === 'rag' ? 'RAG (Strict)' : 'General (Creative)'}
                    </Label>
                    <Switch
                        id="mode-switch"
                        checked={mode === 'general'}
                        onCheckedChange={(c) => setMode(c ? 'general' : 'rag')}
                        className="scale-75 data-[state=checked]:bg-purple-600"
                    />
                    {onClose && (
                        <Button variant="ghost" size="sm" className="h-6 w-6 p-0 ml-2 text-slate-400 hover:text-white" onClick={onClose}>
                            <X className="w-4 h-4" />
                        </Button>
                    )}
                </div>
            </CardHeader>

            <CardContent className="flex-1 p-0 flex flex-col min-h-0">
                <ScrollArea className="flex-1 p-4" ref={scrollRef}>
                    <div className="space-y-4">
                        {messages.length === 0 && (
                            <div className="text-center text-slate-500 py-8 text-xs">
                                <Sparkles className="w-8 h-8 mx-auto mb-2 opacity-20" />
                                <p>Ask about this candidate or general ADC questions.</p>
                                <p className="mt-1 opacity-50">"What is the target mechanism?"</p>
                            </div>
                        )}
                        {messages.map((msg, idx) => (
                            <div key={idx} className={`flex gap-3 ${msg.role === 'user' ? 'flex-row-reverse' : ''}`}>
                                <div className={`w-6 h-6 rounded-full flex items-center justify-center shrink-0 ${msg.role === 'user' ? 'bg-slate-700' : (msg.mode === 'general' ? 'bg-purple-900/50' : 'bg-blue-900/50')
                                    }`}>
                                    {msg.role === 'user' ? <User className="w-3 h-3" /> : <Bot className="w-3 h-3" />}
                                </div>
                                <div className={`flex-1 rounded-lg p-3 text-sm ${msg.role === 'user' ? 'bg-slate-800' : 'bg-slate-950 border border-slate-800'
                                    }`}>
                                    {renderContent(msg)}
                                </div>
                            </div>
                        ))}
                        {isLoading && (
                            <div className="flex gap-3">
                                <div className="w-6 h-6 rounded-full bg-blue-900/50 flex items-center justify-center shrink-0">
                                    <Bot className="w-3 h-3 animate-pulse" />
                                </div>
                                <div className="bg-slate-950 border border-slate-800 rounded-lg p-3">
                                    <span className="text-xs text-slate-500 animate-pulse">Thinking...</span>
                                </div>
                            </div>
                        )}
                    </div>
                </ScrollArea>

                <div className="p-3 border-t border-slate-800 bg-slate-950">
                    <form
                        onSubmit={(e) => {
                            e.preventDefault()
                            handleSend(input)
                        }}
                        className="flex gap-2"
                    >
                        <Input
                            value={input}
                            onChange={e => setInput(e.target.value)}
                            placeholder={mode === 'rag' ? "Ask using internal DB..." : "Ask general questions..."}
                            className="bg-slate-900 border-slate-800 text-xs h-9"
                        />
                        <Button type="submit" size="sm" disabled={isLoading || !input.trim()} className="h-9 w-9 p-0">
                            <ExternalLink className="w-4 h-4" />
                        </Button>
                    </form>
                </div>
            </CardContent>
        </Card>
    )
}
