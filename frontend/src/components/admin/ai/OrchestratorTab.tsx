import { useState } from 'react'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Switch } from '@/components/ui/switch'
import { Label } from '@/components/ui/label'
import {
    Network,
    Cpu,
    Save
} from 'lucide-react'
import { toast } from 'sonner'

export function OrchestratorTab() {
    const [settings, setSettings] = useState({
        autoRouting: true,
        parallelExecution: false
    })

    const handleSave = () => {
        toast.success('Orchestrator settings saved')
    }

    return (
        <div className="grid lg:grid-cols-2 gap-6">
            <Card className="bg-slate-900 border-slate-800">
                <CardHeader>
                    <CardTitle className="text-white flex items-center gap-2">
                        <Network className="w-5 h-5" />
                        Global Routing Settings
                    </CardTitle>
                    <CardDescription className="text-slate-400">
                        Configure how the main agent router delegates tasks.
                    </CardDescription>
                </CardHeader>
                <CardContent className="space-y-6">
                    <div className="flex items-center justify-between">
                        <div className="space-y-0.5">
                            <Label className="text-white">Auto-Routing</Label>
                            <p className="text-sm text-slate-400">Automatically select best agent for query</p>
                        </div>
                        <Switch
                            checked={settings.autoRouting}
                            onCheckedChange={(c) => setSettings(s => ({ ...s, autoRouting: c }))}
                        />
                    </div>
                    <div className="flex items-center justify-between">
                        <div className="space-y-0.5">
                            <Label className="text-white">Parallel Execution</Label>
                            <p className="text-sm text-slate-400">Run multiple agents simultaneously (Higher Cost)</p>
                        </div>
                        <Switch
                            checked={settings.parallelExecution}
                            onCheckedChange={(c) => setSettings(s => ({ ...s, parallelExecution: c }))}
                        />
                    </div>
                </CardContent>
            </Card>

            <Card className="bg-slate-900 border-slate-800">
                <CardHeader>
                    <CardTitle className="text-white flex items-center gap-2">
                        <Cpu className="w-5 h-5" />
                        Model Parameters
                    </CardTitle>
                    <CardDescription className="text-slate-400">
                        Default parameters for all agents.
                    </CardDescription>
                </CardHeader>
                <CardContent className="space-y-6">
                    <div className="space-y-2">
                        <Label className="text-white">Default Model</Label>
                        <div className="p-3 bg-slate-800 rounded-lg text-slate-300 text-sm">
                            GPT-4o (Fixed)
                        </div>
                    </div>
                    <Button onClick={handleSave} className="w-full bg-purple-500 hover:bg-purple-600">
                        <Save className="w-4 h-4 mr-2" />
                        Save Configuration
                    </Button>
                </CardContent>
            </Card>
        </div>
    )
}
