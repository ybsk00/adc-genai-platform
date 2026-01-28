import { useState } from 'react'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import { Button } from '@/components/ui/button'
import { TrendingUp, RefreshCw, ExternalLink } from 'lucide-react'
import { motion, AnimatePresence } from 'framer-motion'

interface NewsItem {
    id: string
    title: string
    source: string
    date: string
    type: 'clinical' | 'business' | 'regulatory'
}

const initialNews: NewsItem[] = [
    {
        id: '1',
        title: 'Pfizer announces new linker technology for next-gen ADCs',
        source: 'Perplexity AI',
        date: '2026-01-17 14:30',
        type: 'business'
    },
    {
        id: '2',
        title: 'Enhertu receives FDA priority review for HER2-low breast cancer',
        source: 'ClinicalTrials.gov',
        date: '2026-01-17 11:00',
        type: 'regulatory'
    },
    {
        id: '3',
        title: 'Phase 1 trial of TROP2-ADC shows promising safety profile',
        source: 'PubMed',
        date: '2026-01-16 09:15',
        type: 'clinical'
    }
]

export function RealTimeFeed() {
    const [news, setNews] = useState<NewsItem[]>(initialNews)
    const [isRefreshing, setIsRefreshing] = useState(false)

    const handleRefresh = async () => {
        setIsRefreshing(true)
        // Simulate API call
        await new Promise(resolve => setTimeout(resolve, 1500))

        // Add a new mock item
        const newItem: NewsItem = {
            id: Math.random().toString(),
            title: 'New competitor analysis report available for lung cancer targets',
            source: 'Market Intelligence',
            date: new Date().toISOString().slice(0, 16).replace('T', ' '),
            type: 'business'
        }

        setNews(prev => [newItem, ...prev.slice(0, 4)])
        setIsRefreshing(false)
    }

    const getTypeColor = (type: string) => {
        switch (type) {
            case 'clinical': return 'bg-blue-500/10 text-blue-400 border-blue-500/20'
            case 'business': return 'bg-purple-500/10 text-purple-400 border-purple-500/20'
            case 'regulatory': return 'bg-green-500/10 text-green-400 border-green-500/20'
            default: return 'bg-slate-800 text-slate-400 border-slate-700'
        }
    }

    return (
        <Card className="bg-slate-900 border-slate-800">
            <CardHeader className="flex flex-row items-center justify-between pb-2">
                <div className="space-y-1">
                    <CardTitle className="flex items-center gap-2 text-base text-white">
                        <TrendingUp className="w-5 h-5 text-blue-500" />
                        Real-time ADC Feed
                    </CardTitle>
                    <CardDescription className="text-slate-400">AI-curated updates from 50+ sources</CardDescription>
                </div>
                <Button
                    variant="ghost"
                    size="icon"
                    onClick={handleRefresh}
                    disabled={isRefreshing}
                    className={`text-slate-400 hover:text-white hover:bg-slate-800 ${isRefreshing ? 'animate-spin' : ''}`}
                >
                    <RefreshCw className="w-4 h-4" />
                </Button>
            </CardHeader>
            <CardContent>
                <div className="space-y-4">
                    <AnimatePresence mode="popLayout">
                        {news.map((item) => (
                            <motion.div
                                key={item.id}
                                initial={{ opacity: 0, y: -20 }}
                                animate={{ opacity: 1, y: 0 }}
                                exit={{ opacity: 0, scale: 0.95 }}
                                layout
                                className="p-3 rounded-lg bg-slate-950 border border-slate-800 hover:border-slate-700 transition-colors"
                            >
                                <div className="flex justify-between items-start gap-2">
                                    <p className="text-sm font-medium text-slate-200 leading-snug">
                                        {item.title}
                                    </p>
                                    <ExternalLink className="w-3 h-3 text-slate-500 flex-shrink-0 mt-1" />
                                </div>
                                <div className="flex items-center gap-2 mt-2">
                                    <Badge variant="outline" className={`text-[10px] px-1.5 py-0 h-5 ${getTypeColor(item.type)}`}>
                                        {item.type}
                                    </Badge>
                                    <span className="text-[10px] text-slate-500">
                                        {item.source} • {item.date}
                                    </span>
                                </div>
                            </motion.div>
                        ))}
                    </AnimatePresence>
                </div>
            </CardContent>
        </Card>
    )
}
