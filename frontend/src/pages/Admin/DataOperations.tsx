import { useState } from 'react'
import { motion } from 'framer-motion'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import {
    Database,
    FileText,
    Table as TableIcon,
    Layers,
    Bot
} from 'lucide-react'
import { DataSourcesTab } from '@/components/admin/data/DataSourcesTab'
import { StagingAreaTab } from '@/components/admin/data/StagingAreaTab'
import { GoldenSetLibraryTab } from '@/components/admin/data/GoldenSetLibraryTab'
import { KnowledgeBaseTab } from '@/components/admin/data/KnowledgeBaseTab'
import { CommercialCrawlerTab } from '@/components/admin/data/CommercialCrawlerTab'

export function DataOperations() {
    return (
        <div className="space-y-6">
            {/* Page Header */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
            >
                <h1 className="text-2xl font-bold text-white">Data Operations</h1>
                <p className="text-slate-400 mt-1">Manage data pipelines, golden sets, and knowledge base.</p>
            </motion.div>

            <Tabs defaultValue="sources" className="space-y-4">
                <TabsList className="bg-slate-900 border border-slate-800">
                    <TabsTrigger value="sources" className="data-[state=active]:bg-slate-800 text-slate-400 data-[state=active]:text-white">
                        <Database className="w-4 h-4 mr-2" />
                        Data Sources
                    </TabsTrigger>
                    <TabsTrigger value="staging" className="data-[state=active]:bg-slate-800 text-slate-400 data-[state=active]:text-white">
                        <Layers className="w-4 h-4 mr-2" />
                        Staging Area
                    </TabsTrigger>
                    <TabsTrigger value="library" className="data-[state=active]:bg-slate-800 text-slate-400 data-[state=active]:text-white">
                        <TableIcon className="w-4 h-4 mr-2" />
                        Golden Set Library
                    </TabsTrigger>
                    <TabsTrigger value="knowledge" className="data-[state=active]:bg-slate-800 text-slate-400 data-[state=active]:text-white">
                        <FileText className="w-4 h-4 mr-2" />
                        Knowledge Base
                    </TabsTrigger>
                    <TabsTrigger value="crawler" className="data-[state=active]:bg-slate-800 text-slate-400 data-[state=active]:text-white">
                        <Bot className="w-4 h-4 mr-2" />
                        Commercial DB
                    </TabsTrigger>
                </TabsList>

                <TabsContent value="sources">
                    <DataSourcesTab />
                </TabsContent>

                <TabsContent value="staging">
                    <StagingAreaTab />
                </TabsContent>

                <TabsContent value="library">
                    <GoldenSetLibraryTab />
                </TabsContent>

                <TabsContent value="knowledge">
                    <KnowledgeBaseTab />
                </TabsContent>

                <TabsContent value="crawler">
                    <CommercialCrawlerTab />
                </TabsContent>
            </Tabs>
        </div>
    )
}

