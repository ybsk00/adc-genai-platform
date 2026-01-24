import { useState } from 'react'
import { motion } from 'framer-motion'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Switch } from '@/components/ui/switch'
import {
    Database,
    FlaskConical,
    BookOpen,
    Filter,
} from 'lucide-react'
import { CommercialTable } from '@/components/admin/data/CommercialTable'
import { KnowledgeTable } from '@/components/admin/data/KnowledgeTable'
import { DataSourcesTab } from '@/components/admin/data/DataSourcesTab'

export default function TotalDataInventory() {
    const [activeTab, setActiveTab] = useState('datasources')
    const [missingDataOnly, setMissingDataOnly] = useState(false)
    
    return (
        <div className="space-y-6 p-6 min-h-screen bg-slate-950 text-slate-100">
            {/* 1. Header */}
            <div className="flex flex-col md:flex-row gap-6 items-start justify-between">
                <div>
                    <h1 className="text-3xl font-bold bg-gradient-to-r from-blue-400 to-purple-400 bg-clip-text text-transparent flex items-center gap-3">
                        <Database className="w-8 h-8 text-blue-400" />
                        Total Data Inventory
                    </h1>
                    <p className="text-slate-400 mt-2">
                        Manage Data Sources, Commercial Reagents, and Knowledge Base.
                    </p>
                </div>
            </div>

            {/* 2. Main Multi-Table View */}
            <Card className="bg-slate-900 border-slate-800">
                <CardHeader className="pb-2">
                    <div className="flex justify-between items-center">
                        <CardTitle>Data Operations Center</CardTitle>
                        <div className="flex gap-2">
                            {/* Global Filters (Only show for inventory tabs) */}
                            {activeTab !== 'datasources' && (
                                <div className="flex items-center gap-2 bg-slate-800 rounded-lg px-3 py-1 border border-slate-700">
                                    <Filter className="w-4 h-4 text-slate-400" />
                                    <span className="text-sm text-slate-300">Missing Data Only</span>
                                    <Switch 
                                        id="missing-filter" 
                                        checked={missingDataOnly}
                                        onCheckedChange={setMissingDataOnly}
                                    />
                                </div>
                            )}
                        </div>
                    </div>
                </CardHeader>
                <CardContent>
                    <Tabs value={activeTab} onValueChange={setActiveTab} className="w-full">
                        <TabsList className="grid w-full grid-cols-3 bg-slate-800 mb-6">
                            <TabsTrigger value="datasources" className="data-[state=active]:bg-blue-900/50 data-[state=active]:text-blue-200">
                                <Database className="w-4 h-4 mr-2" />
                                Data Sources (Crawlers)
                            </TabsTrigger>
                            <TabsTrigger value="commercial" className="data-[state=active]:bg-pink-900/50 data-[state=active]:text-pink-200">
                                <FlaskConical className="w-4 h-4 mr-2" />
                                Commercial Reagents
                            </TabsTrigger>
                            <TabsTrigger value="knowledge" className="data-[state=active]:bg-green-900/50 data-[state=active]:text-green-200">
                                <BookOpen className="w-4 h-4 mr-2" />
                                Knowledge Base
                            </TabsTrigger>
                        </TabsList>

                        <TabsContent value="datasources">
                            <DataSourcesTab />
                        </TabsContent>

                        <TabsContent value="commercial">
                            <CommercialTable missingDataOnly={missingDataOnly} />
                        </TabsContent>
                        
                        <TabsContent value="knowledge">
                            <KnowledgeTable missingDataOnly={missingDataOnly} />
                        </TabsContent>
                    </Tabs>
                </CardContent>
            </Card>
        </div>
    )
}

