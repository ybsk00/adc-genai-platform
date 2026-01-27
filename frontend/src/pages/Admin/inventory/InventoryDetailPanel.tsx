import { useState, useRef, useEffect } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { ScrollArea } from '@/components/ui/scroll-area'
import { X, Send, Bot, User, Sparkles, Loader2 } from 'lucide-react'
import { cn } from '@/lib/utils'


interface InventoryDetailPanelProps {
    item: any
    type: 'antibodies' | 'reagents'
    onClose: () => void
}

export function InventoryDetailPanel({ item, type, onClose }: InventoryDetailPanelProps) {
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

    return (
        <Card className="h-full bg-slate-900 border-slate-800 flex flex-col overflow-hidden">
            <CardHeader className="flex flex-row items-center justify-between py-4 border-b border-slate-800 bg-slate-900/50">
                <CardTitle className="text-lg font-semibold truncate pr-4">
                    {item.name || item.id}
                </CardTitle>
                <Button variant="ghost" size="icon" onClick={onClose} className="h-8 w-8 text-slate-400 hover:text-white">
                    <X className="w-4 h-4" />
                </Button>
            </CardHeader>

            <div className="flex-1 overflow-y-auto p-4 space-y-6">
                {/* Details Section */}
                <div className="space-y-4">
                    <h3 className="text-sm font-medium text-slate-400 uppercase tracking-wider">Details</h3>
                    <div className="grid gap-3 text-sm">
                        {type === 'antibodies' ? (
                            <>
                                <DetailRow label="ID" value={item.id} />
                                <DetailRow label="Target" value={item.target} />
                                <DetailRow label="Host" value={item.host_species} />
                                <DetailRow label="Clonality" value={item.clonality} />
                                <DetailRow label="Isotype" value={item.isotype} />
                                <DetailRow label="Conjugate" value={item.conjugate} />
                            </>
                        ) : (
                            <>
                                <DetailRow label="Name" value={item.name} />
                                <DetailRow label="CAS" value={item.cas_number} />
                                <DetailRow label="Source" value={item.source} />
                                <DetailRow label="Purity" value={item.purity} />
                                <DetailRow label="SMILES" value={item.smiles} className="font-mono text-xs break-all" />
                            </>
                        )}
                    </div>
                </div>

                {/* AI Assistant Section */}
                <div className="pt-4 border-t border-slate-800">
                    <h3 className="text-sm font-medium text-slate-400 uppercase tracking-wider mb-4 flex items-center gap-2">
                        <Sparkles className="w-4 h-4 text-purple-400" />
                        Ask AI Assistant
                    </h3>
                    <AIAssistantPanel item={item} type={type} />
                </div>
            </div>
        </Card>
    )
}

function DetailRow({ label, value, className }: { label: string, value: any, className?: string }) {
    if (!value) return null
    return (
        <div className="grid grid-cols-3 gap-2">
            <span className="text-slate-500">{label}</span>
            <span className={cn("col-span-2 text-slate-200", className)}>{value}</span>
        </div>
    )
}

function AIAssistantPanel({ item, type }: { item: any, type: string }) {
    const [messages, setMessages] = useState<{ role: 'user' | 'assistant', content: string }[]>([])
    const [input, setInput] = useState('')
    const [loading, setLoading] = useState(false)
    const scrollRef = useRef<HTMLDivElement>(null)

    useEffect(() => {
        // Reset chat when item changes
        setMessages([])
    }, [item.id])

    useEffect(() => {
        if (scrollRef.current) {
            scrollRef.current.scrollTop = scrollRef.current.scrollHeight
        }
    }, [messages])

    const handleSend = async () => {
        if (!input.trim() || loading) return

        const userMsg = input
        setMessages(prev => [...prev, { role: 'user', content: userMsg }])
        setInput('')
        setLoading(true)

        try {
            const res = await fetch('/api/admin/ai/chat', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    record_id: item.id,
                    message: userMsg,
                    context: item, // Send full item context
                    mode: 'rag'
                })
            })

            if (res.ok) {
                const data = await res.json()
                setMessages(prev => [...prev, { role: 'assistant', content: data.answer }])
            } else {
                setMessages(prev => [...prev, { role: 'assistant', content: "Sorry, I encountered an error." }])
            }
        } catch (error) {
            console.error("AI Chat Error", error)
            setMessages(prev => [...prev, { role: 'assistant', content: "Network error. Please try again." }])
        } finally {
            setLoading(false)
        }
    }

    return (
        <div className="bg-slate-950 rounded-lg border border-slate-800 flex flex-col h-[400px]">
            <ScrollArea className="flex-1 p-4" ref={scrollRef}>
                <div className="space-y-4">
                    {messages.length === 0 && (
                        <div className="text-center text-slate-500 text-sm py-8">
                            <Bot className="w-8 h-8 mx-auto mb-2 opacity-50" />
                            <p>Ask anything about this {type === 'antibodies' ? 'antibody' : 'reagent'}.</p>
                            <p className="text-xs mt-1">Powered by Gemini 2.5 Flash</p>
                        </div>
                    )}
                    {messages.map((msg, idx) => (
                        <div key={idx} className={cn("flex gap-3", msg.role === 'user' ? "justify-end" : "justify-start")}>
                                <Bot className="w-3 h-3 text-purple-400" />
                            </div>
                            <div className="bg-slate-800 rounded-lg px-3 py-2">
                                <Loader2 className="w-4 h-4 animate-spin text-slate-400" />
                            </div>
                        </div>
                    )}
        </div>
            </ScrollArea >
        <div className="p-3 border-t border-slate-800 flex gap-2">
            <Input
                value={input}
                onChange={(e) => setInput(e.target.value)}
                onKeyDown={(e) => e.key === 'Enter' && handleSend()}
                placeholder="Ask a question..."
                className="bg-slate-900 border-slate-700 text-slate-200 focus:ring-purple-500"
                disabled={loading}
            />
            <Button size="icon" onClick={handleSend} disabled={loading || !input.trim()} className="bg-purple-600 hover:bg-purple-700">
                <Send className="w-4 h-4" />
            </Button>
        </div>
        </div >
    )
}
