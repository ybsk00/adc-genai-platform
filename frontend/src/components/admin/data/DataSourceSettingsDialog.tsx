import { useState } from 'react'
import { Dialog, DialogContent, DialogHeader, DialogTitle, DialogDescription, DialogFooter } from '@/components/ui/dialog'
import { Button } from '@/components/ui/button'
import { Label } from '@/components/ui/label'
import { Input } from '@/components/ui/input'
import { Switch } from '@/components/ui/switch'
import { toast } from 'sonner'

interface DataSourceSettingsDialogProps {
    open: boolean
    onOpenChange: (open: boolean) => void
    sourceId: string
    sourceName: string
}

export function DataSourceSettingsDialog({ open, onOpenChange, sourceId, sourceName }: DataSourceSettingsDialogProps) {
    const [autoSync, setAutoSync] = useState(false)
    const [syncInterval, setSyncInterval] = useState('24')

    const handleSave = () => {
        toast.success(`${sourceName} settings saved.`)
        onOpenChange(false)
    }

    return (
        <Dialog open={open} onOpenChange={onOpenChange}>
            <DialogContent className="bg-slate-900 border-slate-800 text-white">
                <DialogHeader>
                    <DialogTitle>{sourceName} Settings</DialogTitle>
                    <DialogDescription className="text-slate-400">
                        Configure synchronization preferences for this data source.
                    </DialogDescription>
                </DialogHeader>
                <div className="grid gap-4 py-4">
                    <div className="flex items-center justify-between">
                        <Label htmlFor="auto-sync" className="text-white">Auto-sync enabled</Label>
                        <Switch
                            id="auto-sync"
                            checked={autoSync}
                            onCheckedChange={setAutoSync}
                        />
                    </div>
                    <div className="grid gap-2">
                        <Label htmlFor="interval" className="text-white">Sync Interval (hours)</Label>
                        <Input
                            id="interval"
                            value={syncInterval}
                            onChange={(e) => setSyncInterval(e.target.value)}
                            className="bg-slate-800 border-slate-700 text-white"
                        />
                    </div>
                </div>
                <DialogFooter>
                    <Button variant="outline" onClick={() => onOpenChange(false)} className="border-slate-700 text-slate-300">
                        Cancel
                    </Button>
                    <Button onClick={handleSave} className="bg-purple-500 hover:bg-purple-600">
                        Save Changes
                    </Button>
                </DialogFooter>
            </DialogContent>
        </Dialog>
    )
}
