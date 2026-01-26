import { useState } from 'react'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import { Card, CardContent } from '@/components/ui/card'
import { FlaskConical, BookOpen, Microscope } from 'lucide-react'
import { StagingAreaTab } from '@/components/admin/data/StagingAreaTab'
import { KnowledgeEvidenceTab } from '@/components/admin/data/KnowledgeEvidenceTab'

export function StagingAreaLayout() {
    const [activeTab, setActiveTab] = useState("candidates")
    const [globalSearchQuery, setGlobalSearchQuery] = useState("")

    const handleNctClick = (nctId: string) => {
        setGlobalSearchQuery(nctId)
        setActiveTab("candidates")
    }

    return (
        <div className="space-y-6">
            <div>
                <h1 className="text-3xl font-bold text-white flex items-center gap-3">
                    <FlaskConical className="w-8 h-8 text-pink-500" />
                    Staging Area
                </h1>
                <p className="text-slate-400 mt-2">
                    Review and promote AI-refined candidates to the Golden Set Library.
                </p>
            </div>

            <Tabs value={activeTab} onValueChange={setActiveTab} className="w-full">
                <TabsList className="grid w-full grid-cols-2 bg-slate-900 border border-slate-800 mb-6">
                    <TabsTrigger
                        value="candidates"
                        className="data-[state=active]:bg-pink-900/20 data-[state=active]:text-pink-200 data-[state=active]:border-pink-500/50 border-b-2 border-transparent transition-all py-3"
                    >
                        <Microscope className="w-4 h-4 mr-2" />
                        Golden Set Candidates
                    </TabsTrigger>
                    <TabsTrigger
                        value="knowledge"
                        className="data-[state=active]:bg-blue-900/20 data-[state=active]:text-blue-200 data-[state=active]:border-blue-500/50 border-b-2 border-transparent transition-all py-3"
                    >
                        <BookOpen className="w-4 h-4 mr-2" />
                        Knowledge Base Evidence
                    </TabsTrigger>
                </TabsList>

                <TabsContent value="candidates" className="mt-0">
                    <Card className="bg-slate-900 border-slate-800">
                        <CardContent className="p-0">
                            <StagingAreaTab
                                initialSearchQuery={globalSearchQuery}
                                onSearchClear={() => setGlobalSearchQuery("")}
                                onSwitchToKnowledge={() => setActiveTab("knowledge")}
                            />
                        </CardContent>
                    </Card>
                </TabsContent>

                <TabsContent value="knowledge" className="mt-0">
                    <Card className="bg-slate-900 border-slate-800">
                        <CardContent className="p-0">
                            <KnowledgeEvidenceTab onNctClick={handleNctClick} />
                        </CardContent>
                    </Card>
                </TabsContent>
            </Tabs>
        </div>
    )
}
