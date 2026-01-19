import { useState } from 'react'
import { motion } from 'framer-motion'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Textarea } from '@/components/ui/textarea'
import { Input } from '@/components/ui/input'
import { Badge } from '@/components/ui/badge'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import {
    Dna,
    Shield,
    FileText,
    Save,
    Undo,
    Send,
    Loader2,
    History,
    Play
} from 'lucide-react'
import { toast } from 'sonner'

// Agent 프롬프트 타입
interface AgentPrompt {
    id: string
    name: string
    icon: React.ElementType
    currentVersion: string
    draftVersion: string | null
    systemPrompt: string
    draftPrompt: string | null
    lastUpdated: string
}

const initialPrompts: AgentPrompt[] = [
    {
        id: 'structure',
        name: 'Structure Agent',
        icon: Dna,
        currentVersion: 'v1.2',
        draftVersion: null,
        systemPrompt: `You are a Structure Analysis Agent specializing in ADC (Antibody-Drug Conjugate) molecular structures.

Your role:
1. Analyze 3D protein structures from PDB/FASTA input
2. Identify optimal conjugation sites on the antibody
3. Predict structural stability after payload attachment
4. Generate molecular visualization data

Output format:
- Conjugation sites: list of residue positions
- Stability score: 0-100
- Risk factors: list of potential issues`,
        draftPrompt: null,
        lastUpdated: '2026-01-15 14:30'
    },
    {
        id: 'toxicology',
        name: 'Toxicology Agent',
        icon: Shield,
        currentVersion: 'v1.0',
        draftVersion: 'v1.1 (Draft)',
        systemPrompt: `You are a Toxicology Prediction Agent for ADC development.

Your role:
1. Analyze payload toxicity profiles
2. Predict off-target effects based on linker chemistry
3. Compare with known toxicity data from similar compounds
4. Generate risk assessment scores`,
        draftPrompt: `You are a Toxicology Prediction Agent for ADC development.

Your role:
1. Analyze payload toxicity profiles using SMILES strings
2. Predict off-target effects based on linker chemistry and DAR
3. Compare with known toxicity data from similar compounds
4. Generate risk assessment scores with confidence intervals
5. [NEW] Flag specific organ toxicity risks (liver, bone marrow, cardiac)

Output must include:
- Overall toxicity score (0-100)
- Organ-specific risks
- Recommended preclinical tests`,
        lastUpdated: '2026-01-17 09:00'
    },
    {
        id: 'report',
        name: 'Report Writer',
        icon: FileText,
        currentVersion: 'v2.0',
        draftVersion: null,
        systemPrompt: `You are a Report Generator for ADC analysis results.

Your role:
1. Compile findings from all analysis agents
2. Generate executive summary in Korean
3. Create structured sections for each analysis domain
4. Highlight key risks and recommendations

Report sections:
- 개요 (Executive Summary)
- 구조 분석 (Structure Analysis)
- 독성 평가 (Toxicity Assessment)
- 특허 현황 (Patent Landscape)
- 경쟁사 분석 (Competitive Analysis)
- 권고사항 (Recommendations)`,
        draftPrompt: null,
        lastUpdated: '2026-01-10 11:00'
    },
]

// Sandbox 채팅 메시지
interface ChatMessage {
    role: 'user' | 'assistant'
    content: string
}

export function PromptEditorTab() {
    const [prompts, setPrompts] = useState<AgentPrompt[]>(initialPrompts)
    const [activeAgent, setActiveAgent] = useState('structure')
    const [isSaving, setIsSaving] = useState(false)

    // Sandbox state
    const [sandboxInput, setSandboxInput] = useState('')
    const [sandboxMessages, setSandboxMessages] = useState<ChatMessage[]>([])
    const [isGenerating, setIsGenerating] = useState(false)

    const activePrompt = prompts.find(p => p.id === activeAgent)

    const handleEditPrompt = (prompt: string) => {
        setPrompts(prev => prev.map(p =>
            p.id === activeAgent
                ? { ...p, draftPrompt: prompt, draftVersion: `v${parseFloat(p.currentVersion.slice(1)) + 0.1} (Draft)` }
                : p
        ))
    }

    const handleSaveDraft = async () => {
        setIsSaving(true)
        // TODO: API Call
        await new Promise(r => setTimeout(r, 1000))
        toast.success('Draft saved successfully.')
        setIsSaving(false)
    }

    const handlePublish = async () => {
        if (!activePrompt?.draftPrompt) return
        setIsSaving(true)
        // TODO: API Call
        await new Promise(r => setTimeout(r, 1500))

        setPrompts(prev => prev.map(p =>
            p.id === activeAgent
                ? {
                    ...p,
                    systemPrompt: p.draftPrompt || p.systemPrompt,
                    currentVersion: p.draftVersion?.replace(' (Draft)', '') || p.currentVersion,
                    draftPrompt: null,
                    draftVersion: null,
                    lastUpdated: new Date().toISOString().slice(0, 16).replace('T', ' ')
                }
                : p
        ))
        toast.success('Prompt published to Live!')
        setIsSaving(false)
    }

    const handleRevert = () => {
        setPrompts(prev => prev.map(p =>
            p.id === activeAgent
                ? { ...p, draftPrompt: null, draftVersion: null }
                : p
        ))
        toast.info('Changes reverted.')
    }

    const handleSandboxSend = async () => {
        if (!sandboxInput.trim()) return

        const userMessage: ChatMessage = { role: 'user', content: sandboxInput }
        setSandboxMessages(prev => [...prev, userMessage])
        setSandboxInput('')
        setIsGenerating(true)

        // TODO: API Call
        await new Promise(r => setTimeout(r, 2000))

        const assistantMessage: ChatMessage = {
            role: 'assistant',
            content: `[${activePrompt?.name} Response]\n\nAnalysis for LIV-1:\n- Structural Stability: 87/100\n- Recommended DAR: 4\n- Key Risk: Aggregation potential (Medium)\n\nPlease ask if you need detailed residue analysis.`
        }

        setSandboxMessages(prev => [...prev, assistantMessage])
        setIsGenerating(false)
    }

    return (
        <Tabs value={activeAgent} onValueChange={setActiveAgent} className="space-y-4">
            <TabsList className="bg-slate-900 border border-slate-800">
                {prompts.map((prompt) => (
                    <TabsTrigger
                        key={prompt.id}
                        value={prompt.id}
                        className="data-[state=active]:bg-purple-500/20 data-[state=active]:text-white"
                    >
                        <prompt.icon className="w-4 h-4 mr-2" />
                        {prompt.name}
                        {prompt.draftVersion && (
                            <Badge className="ml-2 bg-amber-500/20 text-amber-400 text-xs">Draft</Badge>
                        )}
                    </TabsTrigger>
                ))}
            </TabsList>

            {prompts.map((prompt) => (
                <TabsContent key={prompt.id} value={prompt.id}>
                    <div className="grid lg:grid-cols-2 gap-6">
                        {/* Prompt Editor */}
                        <motion.div
                            initial={{ opacity: 0, x: -20 }}
                            animate={{ opacity: 1, x: 0 }}
                        >
                            <Card className="bg-slate-900 border-slate-800 h-full">
                                <CardHeader>
                                    <div className="flex items-center justify-between">
                                        <div>
                                            <CardTitle className="text-white flex items-center gap-2">
                                                <prompt.icon className="w-5 h-5" />
                                                System Prompt
                                            </CardTitle>
                                            <CardDescription className="text-slate-400 mt-1">
                                                Last Updated: {prompt.lastUpdated}
                                            </CardDescription>
                                        </div>
                                        <div className="flex items-center gap-2">
                                            <Badge variant="outline" className="border-green-500/30 text-green-400">
                                                {prompt.currentVersion} (Live)
                                            </Badge>
                                            {prompt.draftVersion && (
                                                <Badge variant="outline" className="border-amber-500/30 text-amber-400">
                                                    {prompt.draftVersion}
                                                </Badge>
                                            )}
                                        </div>
                                    </div>
                                </CardHeader>
                                <CardContent className="space-y-4">
                                    <Textarea
                                        value={prompt.draftPrompt || prompt.systemPrompt}
                                        onChange={(e) => handleEditPrompt(e.target.value)}
                                        className="bg-slate-800 border-slate-700 text-white font-mono text-sm min-h-[400px]"
                                        placeholder="System prompt..."
                                    />

                                    <div className="flex items-center justify-between">
                                        <Button
                                            variant="outline"
                                            size="sm"
                                            onClick={handleRevert}
                                            disabled={!prompt.draftPrompt}
                                            className="border-slate-700 text-slate-300"
                                        >
                                            <Undo className="w-4 h-4 mr-2" />
                                            Revert
                                        </Button>
                                        <div className="flex gap-2">
                                            <Button
                                                variant="outline"
                                                size="sm"
                                                onClick={handleSaveDraft}
                                                disabled={isSaving || !prompt.draftPrompt}
                                                className="border-slate-700 text-slate-300"
                                            >
                                                {isSaving ? (
                                                    <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                                                ) : (
                                                    <Save className="w-4 h-4 mr-2" />
                                                )}
                                                Save Draft
                                            </Button>
                                            <Button
                                                size="sm"
                                                onClick={handlePublish}
                                                disabled={isSaving || !prompt.draftPrompt}
                                                className="bg-green-500 hover:bg-green-600"
                                            >
                                                <Play className="w-4 h-4 mr-2" />
                                                Publish to Live
                                            </Button>
                                        </div>
                                    </div>

                                    <Button variant="ghost" size="sm" className="w-full text-slate-400 hover:text-white">
                                        <History className="w-4 h-4 mr-2" />
                                        View Version History
                                    </Button>
                                </CardContent>
                            </Card>
                        </motion.div>

                        {/* Sandbox Test Chat */}
                        <motion.div
                            initial={{ opacity: 0, x: 20 }}
                            animate={{ opacity: 1, x: 0 }}
                        >
                            <Card className="bg-slate-900 border-slate-800 h-full flex flex-col">
                                <CardHeader>
                                    <CardTitle className="text-white">Sandbox Test</CardTitle>
                                    <CardDescription className="text-slate-400">
                                        Test your prompt changes before publishing.
                                    </CardDescription>
                                </CardHeader>
                                <CardContent className="flex-1 flex flex-col">
                                    <div className="flex-1 bg-slate-800/50 rounded-lg p-4 mb-4 min-h-[350px] max-h-[350px] overflow-y-auto">
                                        {sandboxMessages.length === 0 ? (
                                            <div className="h-full flex items-center justify-center text-slate-500 text-sm">
                                                Type a message to test the agent...
                                            </div>
                                        ) : (
                                            <div className="space-y-4">
                                                {sandboxMessages.map((msg, i) => (
                                                    <div
                                                        key={i}
                                                        className={`flex ${msg.role === 'user' ? 'justify-end' : 'justify-start'}`}
                                                    >
                                                        <div className={`max-w-[80%] p-3 rounded-lg ${msg.role === 'user'
                                                            ? 'bg-purple-500/20 text-white'
                                                            : 'bg-slate-700 text-slate-200'
                                                            }`}>
                                                            <pre className="whitespace-pre-wrap text-sm font-sans">{msg.content}</pre>
                                                        </div>
                                                    </div>
                                                ))}
                                                {isGenerating && (
                                                    <div className="flex justify-start">
                                                        <div className="bg-slate-700 p-3 rounded-lg">
                                                            <Loader2 className="w-4 h-4 animate-spin text-purple-400" />
                                                        </div>
                                                    </div>
                                                )}
                                            </div>
                                        )}
                                    </div>

                                    <div className="flex gap-2">
                                        <Input
                                            value={sandboxInput}
                                            onChange={(e) => setSandboxInput(e.target.value)}
                                            onKeyDown={(e) => e.key === 'Enter' && !e.shiftKey && handleSandboxSend()}
                                            placeholder="Type a test message..."
                                            className="bg-slate-800 border-slate-700 text-white"
                                        />
                                        <Button
                                            onClick={handleSandboxSend}
                                            disabled={isGenerating || !sandboxInput.trim()}
                                            className="bg-purple-500 hover:bg-purple-600"
                                        >
                                            <Send className="w-4 h-4" />
                                        </Button>
                                    </div>
                                </CardContent>
                            </Card>
                        </motion.div>
                    </div>
                </TabsContent>
            ))}
        </Tabs>
    )
}
