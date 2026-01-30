import { useState, useEffect } from 'react'
import { Dialog, DialogContent, DialogHeader, DialogTitle, DialogDescription, DialogFooter } from '@/components/ui/dialog'
import { Button } from '@/components/ui/button'
import { Label } from '@/components/ui/label'
import { Input } from '@/components/ui/input'
import { Switch } from '@/components/ui/switch'
import { toast } from 'sonner'
import { API_BASE_URL } from '@/lib/api'

interface DataSourceSettingsDialogProps {
    open: boolean
    onOpenChange: (open: boolean) => void
    sourceId: string
    sourceName: string
}

export function DataSourceSettingsDialog({ open, onOpenChange, sourceId, sourceName }: DataSourceSettingsDialogProps) {
    const [autoSync, setAutoSync] = useState(false)
    const [syncInterval, setSyncInterval] = useState('24')
    const [loading, setLoading] = useState(false)

    useEffect(() => {
        if (open && sourceId) {
            fetchSettings()
        }
    }, [open, sourceId])

    const fetchSettings = async () => {
        try {
            const response = await fetch(`${API_BASE_URL}/api/scheduler/settings/${sourceId}`)
            if (response.ok) {
                const data = await response.json()
                setAutoSync(data.auto_sync)
                setSyncInterval(data.sync_interval_hours.toString())
            }
        } catch (error) {
            console.error('Failed to fetch settings:', error)
        }
    }

    const handleSave = async () => {
        setLoading(true)
        try {
            const response = await fetch(`${API_BASE_URL}/api/scheduler/settings`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    source_id: sourceId,
                    auto_sync: autoSync,
                    sync_interval_hours: parseInt(syncInterval) || 24
                })
            })

            if (response.ok) {
                toast.success(`${sourceName} settings saved.`)
                onOpenChange(false)
            } else {
                toast.error('Failed to save settings.')
            }
        } catch (error) {
            toast.error('Error saving settings.')
        } finally {
            setLoading(false)
        }
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
                            disabled={loading}
                        />
                    </div>
                    <div className="grid gap-2">
                        <Label htmlFor="interval" className="text-white">Sync Interval (hours)</Label>
                        <Input
                            id="interval"
                            type="number"
                            value={syncInterval}
                            onChange={(e) => setSyncInterval(e.target.value)}
                            disabled={loading}
                            className="bg-slate-800 border-slate-700 text-white"
                        />
                    </div>
                </div>
                <DialogFooter>
                    <Button variant="outline" onClick={() => onOpenChange(false)} disabled={loading} className="border-slate-700 text-slate-300">
                        Cancel
                    </Button>
                    <Button onClick={handleSave} disabled={loading} className="bg-purple-500 hover:bg-purple-600">
                        {loading ? 'Saving...' : 'Save Changes'}
                    </Button>
                </DialogFooter>
            </DialogContent>
        </Dialog>
    )
}
