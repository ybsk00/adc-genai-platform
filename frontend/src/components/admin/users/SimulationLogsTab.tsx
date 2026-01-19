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
import { Loader2, RefreshCw, FileText, CheckCircle, XCircle, Clock, FileDown, Eye } from 'lucide-react'
import { toast } from 'sonner'
import { format } from 'date-fns'
import { SimulationDetailDrawer } from './SimulationDetailDrawer'
import { API_BASE_URL } from '@/lib/api'

interface SimulationLog {
    id: string
    title: string
    status: 'completed' | 'failed' | 'processing' | 'pending'
    created_at: string
    report_url: string | null
    input_params: {
        target?: string
        antibody?: string
        payload?: string
        linker?: string
        dar?: number
    } | null
    result_summary: {
        success_score?: number
        toxicity_warnings?: string[]
        binding_affinity?: number
        dar_prediction?: number
        grade?: string
        recommendation?: string
    } | null
    profiles: {
        email: string
    } | null
}

export function SimulationLogsTab() {
    const [logs, setLogs] = useState<SimulationLog[]>([])
    const [loading, setLoading] = useState(true)
    const [selectedSimulation, setSelectedSimulation] = useState<SimulationLog | null>(null)
    const [drawerOpen, setDrawerOpen] = useState(false)

    const fetchLogs = async () => {
        setLoading(true)
        try {
            const response = await fetch(`${API_BASE_URL}/api/admin/simulations`)
            if (!response.ok) throw new Error('Failed to fetch logs')
            const result = await response.json()
            setLogs(result)
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

    const handleViewDetails = (log: SimulationLog) => {
        setSelectedSimulation(log)
        setDrawerOpen(true)
    }

    const handleDownloadPDF = (reportUrl: string) => {
        window.open(reportUrl, '_blank')
    }

    const getTargetPayloadBadge = (log: SimulationLog) => {
        const target = log.input_params?.target || 'N/A'
        const payload = log.input_params?.payload || 'N/A'
        return `${target} + ${payload}`
    }

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
                            <TableHead className="text-slate-400">User</TableHead>
                            <TableHead className="text-slate-400">Project Name</TableHead>
                            <TableHead className="text-slate-400">Target / Payload</TableHead>
                            <TableHead className="text-slate-400">Status</TableHead>
                            <TableHead className="text-slate-400">Date</TableHead>
                            <TableHead className="text-slate-400">Report</TableHead>
                            <TableHead className="text-slate-400 text-right">Actions</TableHead>
                        </TableRow>
                    </TableHeader>
                    <TableBody>
                        {loading ? (
                            <TableRow>
                                <TableCell colSpan={7} className="h-24 text-center">
                                    <div className="flex justify-center items-center text-slate-500">
                                        <Loader2 className="w-6 h-6 animate-spin mr-2" />
                                        Loading logs...
                                    </div>
                                </TableCell>
                            </TableRow>
                        ) : logs.length > 0 ? (
                            logs.map((log) => (
                                <TableRow
                                    key={log.id}
                                    className="border-slate-800 hover:bg-slate-800/50 cursor-pointer"
                                    onClick={() => handleViewDetails(log)}
                                >
                                    <TableCell>
                                        <span className="text-slate-400">{log.profiles?.email || 'Unknown'}</span>
                                    </TableCell>
                                    <TableCell>
                                        <div className="flex items-center gap-2">
                                            <FileText className="w-4 h-4 text-slate-500" />
                                            <span className="font-medium text-white">{log.title}</span>
                                        </div>
                                    </TableCell>
                                    <TableCell>
                                        <Badge variant="outline" className="border-purple-500/30 text-purple-400">
                                            {getTargetPayloadBadge(log)}
                                        </Badge>
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
                                    <TableCell>
                                        {log.report_url ? (
                                            <Button
                                                variant="ghost"
                                                size="sm"
                                                onClick={(e) => {
                                                    e.stopPropagation()
                                                    handleDownloadPDF(log.report_url!)
                                                }}
                                                className="text-blue-400 hover:text-blue-300"
                                            >
                                                <FileDown className="w-4 h-4" />
                                            </Button>
                                        ) : (
                                            <span className="text-slate-600 text-sm">-</span>
                                        )}
                                    </TableCell>
                                    <TableCell className="text-right">
                                        <Button
                                            variant="ghost"
                                            size="sm"
                                            onClick={(e) => {
                                                e.stopPropagation()
                                                handleViewDetails(log)
                                            }}
                                            className="text-slate-400 hover:text-white"
                                        >
                                            <Eye className="w-4 h-4 mr-1" />
                                            View Details
                                        </Button>
                                    </TableCell>
                                </TableRow>
                            ))
                        ) : (
                            <TableRow>
                                <TableCell colSpan={7} className="h-24 text-center text-slate-500">
                                    No logs found.
                                </TableCell>
                            </TableRow>
                        )}
                    </TableBody>
                </Table>
            </div>

            <SimulationDetailDrawer
                simulation={selectedSimulation}
                open={drawerOpen}
                onClose={() => setDrawerOpen(false)}
            />
        </div>
    )
}
