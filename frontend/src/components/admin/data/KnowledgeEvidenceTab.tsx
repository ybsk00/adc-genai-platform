import { useState } from 'react'
import { useQuery } from '@tanstack/react-query'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Input } from '@/components/ui/input'
import { Badge } from '@/components/ui/badge'
import { Button } from '@/components/ui/button'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Search, BookOpen, ExternalLink, Loader2, FileText } from 'lucide-react'
import { getSession } from '@/lib/supabase'
import { API_BASE_URL } from '@/lib/api'

interface KnowledgeItem {
    id: number
    source_type: string
    title: string
    summary: string
    relevance_score: number
    rag_status: string
    created_at: string
}

interface KnowledgeEvidenceTabProps {
    onNctClick: (nctId: string) => void
}

export function KnowledgeEvidenceTab({ onNctClick }: KnowledgeEvidenceTabProps) {
    const [searchQuery, setSearchQuery] = useState('')

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

    // Helper to render text with clickable NCT IDs
    const renderContentWithLinks = (text: string) => {
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

    return (
        <div className="h-[calc(100vh-250px)] flex flex-col bg-slate-900 border border-slate-800 rounded-lg overflow-hidden">
            <div className="p-4 border-b border-slate-800 bg-slate-950 flex justify-between items-center">
                <div className="relative w-1/3">
                    <Search className="absolute left-2 top-2.5 h-4 w-4 text-slate-500" />
                    <Input
                        placeholder="Search knowledge base..."
                        className="pl-8 bg-slate-800 border-slate-700 text-white h-9 text-sm"
                        value={searchQuery}
                        onChange={e => setSearchQuery(e.target.value)}
                    />
                </div>
                <div className="text-xs text-slate-500">
                    Showing {filteredItems?.length || 0} items
                </div>
            </div>

            <ScrollArea className="flex-1 p-4">
                {isLoading ? (
                    <div className="flex justify-center p-8"><Loader2 className="animate-spin text-blue-500" /></div>
                ) : (
                    <div className="grid grid-cols-1 gap-4">
                        {filteredItems?.map(item => (
                            <Card key={item.id} className="bg-slate-950 border-slate-800 hover:border-slate-700 transition-colors">
                                <CardHeader className="pb-2">
                                    <div className="flex justify-between items-start">
                                        <div className="flex items-center gap-2">
                                            <Badge variant="outline" className="bg-blue-900/20 text-blue-400 border-blue-800">
                                                {item.source_type}
                                            </Badge>
                                            <span className="text-xs text-slate-500">
                                                {new Date(item.created_at).toLocaleDateString()}
                                            </span>
                                        </div>
                                        <Badge variant="secondary" className={
                                            item.rag_status === 'indexed' ? 'bg-green-900/20 text-green-400' : 'bg-slate-800 text-slate-400'
                                        }>
                                            {item.rag_status}
                                        </Badge>
                                    </div>
                                    <CardTitle className="text-base text-slate-200 mt-2 leading-snug">
                                        {renderContentWithLinks(item.title)}
                                    </CardTitle>
                                </CardHeader>
                                <CardContent>
                                    <p className="text-sm text-slate-400 line-clamp-3">
                                        {renderContentWithLinks(item.summary)}
                                    </p>
                                </CardContent>
                            </Card>
                        ))}
                        {filteredItems?.length === 0 && (
                            <div className="text-center py-12 text-slate-500">
                                <BookOpen className="w-12 h-12 mx-auto mb-4 opacity-20" />
                                <p>No items found matching your search.</p>
                            </div>
                        )}
                    </div>
                )}
            </ScrollArea>
        </div>
    )
}
