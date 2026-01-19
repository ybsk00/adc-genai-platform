import { useState } from 'react'
import { motion } from 'framer-motion'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select'
import { Badge } from '@/components/ui/badge'
import {
    Bot,
    Play,
    ShieldAlert,
    Clock,
    Loader2
} from 'lucide-react'
import { toast } from 'sonner'

export function CommercialCrawlerTab() {
    const [category, setCategory] = useState<string>('Payload')
    const [isRunning, setIsRunning] = useState(false)
    const [logs, setLogs] = useState<string[]>([])

    const handleRunCrawler = async () => {
        setIsRunning(true)
        setLogs(prev => [`[${new Date().toLocaleTimeString()}] Starting crawler for ${category}...`, ...prev])

        try {
            const response = await fetch('/api/admin/crawler/ambeed/run', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ category, max_pages: 1 })
            })

            if (!response.ok) throw new Error('Failed to start crawler')

            const result = await response.json()
            setLogs(prev => [
                `[${new Date().toLocaleTimeString()}] ${result.message}`,
                `[${new Date().toLocaleTimeString()}] Mode: ${result.mode}`,
                ...prev
            ])
            toast.success('Crawler started in background')
        } catch (error) {
            console.error(error)
            toast.error('Failed to start crawler')
            setLogs(prev => [`[${new Date().toLocaleTimeString()}] Error: Failed to start crawler`, ...prev])
        } finally {
            setTimeout(() => setIsRunning(false), 2000) // Simulate "starting" delay
        }
    }

    return (
        <div className="grid lg:grid-cols-3 gap-6">
            {/* Left Column: Control Panel */}
            <div className="lg:col-span-1 space-y-6">
                <Card className="bg-slate-900 border-slate-800">
                    <CardHeader>
                        <CardTitle className="text-white flex items-center gap-2">
                            <Bot className="w-5 h-5 text-purple-400" />
                            Ambeed Stealth Crawler
                        </CardTitle>
                        <CardDescription className="text-slate-400">
                            Securely crawl commercial reagents with anti-blocking technology.
                        </CardDescription>
                    </CardHeader>
                    <CardContent className="space-y-6">
                        <div className="space-y-2">
                            <label className="text-sm text-slate-400">Target Category</label>
                            <Select value={category} onValueChange={setCategory}>
                                <SelectTrigger className="bg-slate-800 border-slate-700 text-white">
                                    <SelectValue />
                                </SelectTrigger>
                                <SelectContent className="bg-slate-800 border-slate-700 text-white">
                                    <SelectItem value="Payload">Payload (Toxins)</SelectItem>
                                    <SelectItem value="Linker">Linkers</SelectItem>
                                    <SelectItem value="Conjugate">Antibody-Drug Conjugates</SelectItem>
                                </SelectContent>
                            </Select>
                        </div>

                        <div className="p-4 bg-slate-800/50 rounded-lg border border-slate-700/50 space-y-3">
                            <div className="flex items-center gap-2 text-sm text-slate-300">
                                <ShieldAlert className="w-4 h-4 text-green-400" />
                                <span>Anti-Blocking Active</span>
                            </div>
                            <div className="flex items-center gap-2 text-sm text-slate-300">
                                <Clock className="w-4 h-4 text-blue-400" />
                                <span>Random Delay (2-5s)</span>
                            </div>
                            <div className="flex items-center gap-2 text-sm text-slate-300">
                                <Bot className="w-4 h-4 text-orange-400" />
                                <span>User-Agent Rotation</span>
                            </div>
                        </div>

                        <Button
                            className="w-full bg-purple-600 hover:bg-purple-700"
                            onClick={handleRunCrawler}
                            disabled={isRunning}
                        >
                            {isRunning ? (
                                <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                            ) : (
                                <Play className="w-4 h-4 mr-2" />
                            )}
                            {isRunning ? 'Starting...' : 'Start Crawling'}
                        </Button>
                    </CardContent>
                </Card>
            </div>

            {/* Right Column: Logs & Status */}
            <div className="lg:col-span-2 space-y-6">
                <Card className="bg-slate-900 border-slate-800 h-full">
                    <CardHeader>
                        <CardTitle className="text-white flex items-center justify-between">
                            <span>Crawler Logs</span>
                            <Badge variant="outline" className="border-slate-700 text-slate-400">
                                Real-time
                            </Badge>
                        </CardTitle>
                    </CardHeader>
                    <CardContent>
                        <div className="h-[400px] bg-black/40 rounded-lg p-4 font-mono text-sm overflow-y-auto border border-slate-800">
                            {logs.length > 0 ? (
                                logs.map((log, i) => (
                                    <div key={i} className="mb-2 text-slate-300 border-b border-slate-800/50 pb-1 last:border-0">
                                        {log}
                                    </div>
                                ))
                            ) : (
                                <div className="h-full flex flex-col items-center justify-center text-slate-500 gap-2">
                                    <Bot className="w-8 h-8 opacity-20" />
                                    <p>Ready to crawl. Select a category and start.</p>
                                </div>
                            )}
                        </div>
                    </CardContent>
                </Card>
            </div>
        </div>
    )
}
