import { useState } from 'react'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import { Card } from '@/components/ui/card'
import { Database, FlaskConical } from 'lucide-react'
import { AntibodyListTab } from './AntibodyListTab'
import { ReagentListTab } from './ReagentListTab'
import { InventoryDetailPanel } from './InventoryDetailPanel'

export default function TotalInventoryLayout() {
    const [activeTab, setActiveTab] = useState('antibodies')
    const [selectedItem, setSelectedItem] = useState<any>(null)

    return (
        <div className="h-[calc(100vh-4rem)] p-6 flex gap-6 bg-slate-950 text-slate-100 overflow-hidden">
            {/* Left Panel: List View (33%) */}
            <div className="w-1/3 flex flex-col min-w-0">
                <div className="mb-6">
                    <h1 className="text-3xl font-bold bg-gradient-to-r from-blue-400 to-purple-400 bg-clip-text text-transparent flex items-center gap-3">
                        <Database className="w-8 h-8 text-blue-400" />
                        ADC Total Inventory
                    </h1>
                    <p className="text-slate-400 mt-2">
                        Manage Antibodies and Reagents with AI-powered analysis.
                    </p>
                </div>

                <Card className="flex-1 bg-slate-900 border-slate-800 flex flex-col min-h-0">
                    <Tabs value={activeTab} onValueChange={(val) => {
                        setActiveTab(val)
                        setSelectedItem(null) // Reset selection on tab change
                    }} className="flex-1 flex flex-col min-h-0">
                        <div className="px-6 pt-6 pb-2">
                            <TabsList className="grid w-full grid-cols-2 bg-slate-800">
                                <TabsTrigger
                                    value="antibodies"
                                    className="data-[state=active]:bg-blue-900/50 data-[state=active]:text-blue-200"
                                >
                                    <Database className="w-4 h-4 mr-2" />
                                    Antibodies
                                </TabsTrigger>
                                <TabsTrigger
                                    value="reagents"
                                    className="data-[state=active]:bg-pink-900/50 data-[state=active]:text-pink-200"
                                >
                                    <FlaskConical className="w-4 h-4 mr-2" />
                                    Reagents
                                </TabsTrigger>
                            </TabsList>
                        </div>

                        <div className="flex-1 min-h-0 overflow-hidden">
                            <TabsContent value="antibodies" className="h-full m-0 border-0">
                                <AntibodyListTab onSelect={setSelectedItem} selectedId={selectedItem?.id} />
                            </TabsContent>
                            <TabsContent value="reagents" className="h-full m-0 border-0">
                                <ReagentListTab onSelect={setSelectedItem} selectedId={selectedItem?.id} />
                            </TabsContent>
                        </div>
                    </Tabs>
                </Card>
            </div>

            {/* Right Panel: Detail View (Rest) */}
            <div className="flex-1 flex flex-col min-w-0">
                <InventoryDetailPanel
                    item={selectedItem}
                    type={activeTab as 'antibodies' | 'reagents'}
                    onClose={() => setSelectedItem(null)}
                />
            </div>
        </div>
    )
}
