import { useState, useEffect } from 'react'
import {
    Table,
    TableBody,
    TableCell,
    TableHead,
    TableHeader,
    TableRow,
} from '@/components/ui/table'
import { Badge } from '@/components/ui/badge'
import { Button } from '@/components/ui/button'
import { Loader2, RefreshCw, FileText, CheckCircle, XCircle, Clock } from 'lucide-react'
import { toast } from 'sonner'
import { format } from 'date-fns'

interface SimulationLog {
    id: string
    name: string
    status: 'completed' | 'failed' | 'processing' | 'pending'
    created_at: string
    profiles: {
        email: string
    } | null
    // Add other fields as needed
}

export function SimulationLogsTab() {
    const [logs, setLogs] = useState<SimulationLog[]>([])
    const [loading, setLoading] = useState(true)

    const fetchLogs = async () => {
        setLoading(true)
        try {
            // Mock Data for now as API might not be fully ready or populated
            // const response = await fetch('/api/admin/simulations')
            // const result = await response.json()

            await new Promise(r => setTimeout(r, 1000))
            const mockLogs: SimulationLog[] = [
                {
                    id: 'sim_1',
                    name: 'HER2-DXd Analysis',
                    status: 'completed',
                    created_at: '2026-01-18T10:00:00Z',
                    profiles: { email: 'kim@biolab.kr' }
                },
                {
                    id: 'sim_2',
                    name: 'TROP2 Linker Optimization',
                    status: 'processing',
                    created_at: '2026-01-18T11:30:00Z',
                    profiles: { email: 'lee@university.edu' }
                },
                {
                    id: 'sim_3',
                    name: 'CD30 Payload Screen',
                    status: 'failed',
                    created_at: '2026-01-17T15:00:00Z',
                    profiles: { email: 'park@pharma.com' }
                }
            ]
            setLogs(mockLogs)
        } catch (error) {
            console.error(error)
            toast.error('Failed to load simulation logs')
        } finally {
            setLoading(false)
        }
    }

    useEffect(() => {
        fetchLogs()
    }, [])

    return (
        <div className="space-y-4">
            <div className="flex justify-end">
                <Button
                    variant="outline"
                    onClick={fetchLogs}
                    disabled={loading}
                    className="border-slate-700 text-slate-300"
                >
                    <RefreshCw className={`w-4 h-4 mr-2 ${loading ? 'animate-spin' : ''}`} />
                    Refresh Logs
                </Button>
            </div>

            <div className="rounded-md border border-slate-800 bg-slate-900">
                <Table>
                    <TableHeader>
                        <TableRow className="border-slate-800 hover:bg-transparent">
                            <TableHead className="text-slate-400">Simulation Name</TableHead>
                            <TableHead className="text-slate-400">User</TableHead>
                            <TableHead className="text-slate-400">Status</TableHead>
                            <TableHead className="text-slate-400">Date</TableHead>
                            <TableHead className="text-slate-400 text-right">Actions</TableHead>
                        </TableRow>
                    </TableHeader>
                    <TableBody>
                        {loading ? (
                            <TableRow>
                                <TableCell colSpan={5} className="h-24 text-center">
                                    <div className="flex justify-center items-center text-slate-500">
                                        <Loader2 className="w-6 h-6 animate-spin mr-2" />
                                        Loading logs...
                                    </div>
                                </TableCell>
                            </TableRow>
                        ) : logs.length > 0 ? (
                            logs.map((log) => (
                                <TableRow key={log.id} className="border-slate-800 hover:bg-slate-800/50">
                                    <TableCell>
                                        <div className="flex items-center gap-2">
                                            <FileText className="w-4 h-4 text-slate-500" />
                                            <span className="font-medium text-white">{log.name}</span>
                                        </div>
                                    </TableCell>
                                    <TableCell>
                                        <span className="text-slate-400">{log.profiles?.email || 'Unknown'}</span>
                                    </TableCell>
                                    <TableCell>
                                        <Badge
                                            variant="outline"
                                            className={
                                                log.status === 'completed' ? 'border-green-500/30 text-green-400' :
                                                    log.status === 'processing' ? 'border-blue-500/30 text-blue-400' :
                                                        log.status === 'failed' ? 'border-red-500/30 text-red-400' :
                                                            'border-slate-600 text-slate-400'
                                            }
                                        >
                                            {log.status === 'completed' && <CheckCircle className="w-3 h-3 mr-1" />}
                                            {log.status === 'processing' && <Loader2 className="w-3 h-3 mr-1 animate-spin" />}
                                            {log.status === 'failed' && <XCircle className="w-3 h-3 mr-1" />}
                                            {log.status}
                                        </Badge>
                                    </TableCell>
                                    <TableCell>
                                        <div className="flex items-center gap-1 text-slate-500 text-sm">
                                            <Clock className="w-3 h-3" />
                                            {format(new Date(log.created_at), 'MMM d, HH:mm')}
                                        </div>
                                    </TableCell>
                                    <TableCell className="text-right">
                                        <Button variant="ghost" size="sm" className="text-slate-400 hover:text-white">
                                            View Details
                                        </Button>
                                    </TableCell>
                                </TableRow>
                            ))
                        ) : (
                            <TableRow>
                                <TableCell colSpan={5} className="h-24 text-center text-slate-500">
                                    No logs found.
                                </TableCell>
                            </TableRow>
                        )}
                    </TableBody>
                </Table>
            </div>
        </div>
    )
}
