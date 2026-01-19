import { useState } from 'react'
import { motion } from 'framer-motion'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Badge } from '@/components/ui/badge'
import { Slider } from '@/components/ui/slider'
import {
    Search,
    Database,
    FileText,
    Settings,
    Loader2
} from 'lucide-react'
import { toast } from 'sonner'

interface RetrievedChunk {
    id: string
    content: string
    source: string
    score: number
}

export function RAGGeneratorTab() {
    const [query, setQuery] = useState('')
    const [isSearching, setIsSearching] = useState(false)
    const [results, setResults] = useState<RetrievedChunk[]>([])
    const [generatedAnswer, setGeneratedAnswer] = useState<string | null>(null)

    // Settings
    const [topK, setTopK] = useState([5])
    const [threshold, setThreshold] = useState([0.7])

    const handleSearch = async () => {
        if (!query.trim()) return
        setIsSearching(true)
        setResults([])
        setGeneratedAnswer(null)

        try {
            // TODO: API Call
            await new Promise(r => setTimeout(r, 1500))

            const mockResults: RetrievedChunk[] = [
                {
                    id: 'c1',
                    content: 'T-DXd demonstrated a confirmed objective response rate of 78.5% in patients with HER2-positive metastatic breast cancer...',
                    source: 'PubMed: Clinical Efficacy of T-DXd (2025)',
                    score: 0.92
                },
                {
                    id: 'c2',
                    content: 'Interstitial lung disease (ILD) remains a significant risk factor for T-DXd, requiring careful monitoring...',
                    source: 'Safety Report: ADC Toxicity Profile',
                    score: 0.88
                },
                {
                    id: 'c3',
                    content: 'The market for HER2-targeted ADCs is expected to reach $15B by 2030...',
                    source: 'Market Analysis Report Q4 2025',
                    score: 0.75
                }
            ]
            setResults(mockResults)
            setGeneratedAnswer("Based on the retrieved documents, T-DXd shows high efficacy (78.5% ORR) in HER2+ breast cancer but carries a risk of ILD. The market outlook remains strong.")
        } catch (error) {
            toast.error('Retrieval failed')
        } finally {
            setIsSearching(false)
        }
    }

    return (
        <div className="grid lg:grid-cols-3 gap-6">
            {/* Left Column: Search & Settings */}
            <div className="lg:col-span-1 space-y-6">
                <Card className="bg-slate-900 border-slate-800">
                    <CardHeader>
                        <CardTitle className="text-white flex items-center gap-2">
                            <Settings className="w-5 h-5" />
                            Retrieval Settings
                        </CardTitle>
                    </CardHeader>
                    <CardContent className="space-y-6">
                        <div className="space-y-2">
                            <div className="flex justify-between text-sm">
                                <span className="text-slate-400">Top K (Chunks)</span>
                                <span className="text-white font-medium">{topK[0]}</span>
                            </div>
                            <Slider
                                value={topK}
                                onValueChange={setTopK}
                                max={20}
                                min={1}
                                step={1}
                                className="py-2"
                            />
                        </div>
                        <div className="space-y-2">
                            <div className="flex justify-between text-sm">
                                <span className="text-slate-400">Similarity Threshold</span>
                                <span className="text-white font-medium">{threshold[0]}</span>
                            </div>
                            <Slider
                                value={threshold}
                                onValueChange={setThreshold}
                                max={1}
                                min={0}
                                step={0.05}
                                className="py-2"
                            />
                        </div>
                    </CardContent>
                </Card>

                <Card className="bg-slate-900 border-slate-800">
                    <CardHeader>
                        <CardTitle className="text-white">Query Input</CardTitle>
                    </CardHeader>
                    <CardContent className="space-y-4">
                        <Input
                            placeholder="Enter your question..."
                            value={query}
                            onChange={(e) => setQuery(e.target.value)}
                            onKeyDown={(e) => e.key === 'Enter' && handleSearch()}
                            className="bg-slate-800 border-slate-700 text-white"
                        />
                        <Button
                            className="w-full bg-purple-500 hover:bg-purple-600"
                            onClick={handleSearch}
                            disabled={isSearching || !query.trim()}
                        >
                            {isSearching ? (
                                <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                            ) : (
                                <Search className="w-4 h-4 mr-2" />
                            )}
                            Retrieve & Generate
                        </Button>
                    </CardContent>
                </Card>
            </div>

            {/* Right Column: Results */}
            <div className="lg:col-span-2 space-y-6">
                {/* Generated Answer */}
                {generatedAnswer && (
                    <motion.div
                        initial={{ opacity: 0, y: 10 }}
                        animate={{ opacity: 1, y: 0 }}
                    >
                        <Card className="bg-slate-900 border-slate-800 border-l-4 border-l-purple-500">
                            <CardHeader>
                                <CardTitle className="text-white flex items-center gap-2">
                                    <Database className="w-5 h-5 text-purple-400" />
                                    Generated Answer
                                </CardTitle>
                            </CardHeader>
                            <CardContent>
                                <p className="text-slate-300 leading-relaxed">
                                    {generatedAnswer}
                                </p>
                            </CardContent>
                        </Card>
                    </motion.div>
                )}

                {/* Retrieved Chunks */}
                <div className="space-y-4">
                    <h3 className="text-lg font-semibold text-white flex items-center gap-2">
                        <FileText className="w-5 h-5" />
                        Retrieved Context
                        <Badge variant="secondary" className="bg-slate-800 text-slate-400">
                            {results.length} chunks
                        </Badge>
                    </h3>

                    {results.length > 0 ? (
                        results.map((chunk, index) => (
                            <motion.div
                                key={chunk.id}
                                initial={{ opacity: 0, x: 10 }}
                                animate={{ opacity: 1, x: 0 }}
                                transition={{ delay: index * 0.1 }}
                            >
                                <Card className="bg-slate-900 border-slate-800 hover:border-slate-700 transition-colors">
                                    <CardContent className="p-4">
                                        <div className="flex justify-between items-start mb-2">
                                            <Badge variant="outline" className="border-blue-500/30 text-blue-400">
                                                Score: {chunk.score}
                                            </Badge>
                                            <span className="text-xs text-slate-500 truncate max-w-[200px]">
                                                {chunk.source}
                                            </span>
                                        </div>
                                        <p className="text-sm text-slate-300">
                                            {chunk.content}
                                        </p>
                                    </CardContent>
                                </Card>
                            </motion.div>
                        ))
                    ) : (
                        <div className="text-center py-12 text-slate-500 border border-dashed border-slate-800 rounded-lg">
                            No context retrieved yet. Run a search to see RAG results.
                        </div>
                    )}
                </div>
            </div>
        </div>
    )
}
