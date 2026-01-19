import { motion } from 'framer-motion'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import {
    Network,
    Bot,
    Database
} from 'lucide-react'
import { OrchestratorTab } from '@/components/admin/ai/OrchestratorTab'
import { PromptEditorTab } from '@/components/admin/ai/PromptEditorTab'
import { RAGGeneratorTab } from '@/components/admin/ai/RAGGeneratorTab'

export function AITuning() {
    return (
        <div className="space-y-6">
            {/* Page Header */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
            >
                <h1 className="text-2xl font-bold text-white">AI Tuning</h1>
                <p className="text-slate-400 mt-1">Manage AI agents, prompts, and RAG pipeline.</p>
            </motion.div>

            <Tabs defaultValue="agents" className="space-y-4">
                <TabsList className="bg-slate-900 border border-slate-800">
                    <TabsTrigger value="orchestrator" className="data-[state=active]:bg-slate-800 text-slate-400 data-[state=active]:text-white">
                        <Network className="w-4 h-4 mr-2" />
                        Orchestrator
                    </TabsTrigger>
                    <TabsTrigger value="agents" className="data-[state=active]:bg-slate-800 text-slate-400 data-[state=active]:text-white">
                        <Bot className="w-4 h-4 mr-2" />
                        Agents & Prompts
                    </TabsTrigger>
                    <TabsTrigger value="rag" className="data-[state=active]:bg-slate-800 text-slate-400 data-[state=active]:text-white">
                        <Database className="w-4 h-4 mr-2" />
                        RAG Generator
                    </TabsTrigger>
                </TabsList>

                <TabsContent value="orchestrator">
                    <OrchestratorTab />
                </TabsContent>

                <TabsContent value="agents">
                    <PromptEditorTab />
                </TabsContent>

                <TabsContent value="rag">
                    <RAGGeneratorTab />
                </TabsContent>
            </Tabs>
        </div>
    )
}
