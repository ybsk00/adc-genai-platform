import { useState, useEffect, useCallback } from 'react'
import { motion } from 'framer-motion'
import { Card, CardContent } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import {
    Activity,
    Loader2
} from 'lucide-react'
import { API_BASE_URL } from '@/lib/api'

// System services to monitor
const systemServices = [
    { id: 'api', name: 'API Server', endpoint: `${API_BASE_URL}/health` },
    { id: 'worker', name: 'Background Worker', endpoint: `${API_BASE_URL}/api/worker/health` },
    { id: 'db', name: 'Supabase DB', endpoint: `${API_BASE_URL}/api/db/health` },
]

type ServiceStatus = 'operational' | 'degraded' | 'down' | 'checking'

function useHealthCheck() {
    const [statuses, setStatuses] = useState<Record<string, ServiceStatus>>(
        Object.fromEntries(systemServices.map(s => [s.id, 'checking']))
    )

    const checkHealth = useCallback(async () => {
        for (const service of systemServices) {
            try {
                // Mocking health check for now as endpoints might not be fully ready
                // In production, uncomment fetch
                // const response = await fetch(service.endpoint)
                // const isHealthy = response.ok
                const isHealthy = true // Mock

                setStatuses(prev => ({
                    ...prev,
                    [service.id]: isHealthy ? 'operational' : 'down'
                }))
            } catch {
                setStatuses(prev => ({
                    ...prev,
                    [service.id]: 'down'
                }))
            }
        }
    }, [])

    useEffect(() => {
        checkHealth()
        const interval = setInterval(checkHealth, 60000)
        return () => clearInterval(interval)
    }, [checkHealth])

    return statuses
}

export function AdminOverview() {
    const serviceStatuses = useHealthCheck()

    return (
        <div className="space-y-6">
            {/* Page Header */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
            >
                <h1 className="text-2xl font-bold text-white">Admin Overview</h1>
                <p className="text-slate-400 mt-1">System Health & Real-time Metrics</p>
            </motion.div>

            {/* System Status - Traffic Lights */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.1 }}
                className="grid grid-cols-1 md:grid-cols-3 gap-4"
            >
                {systemServices.map((service) => {
                    const status = serviceStatuses[service.id]
                    return (
                        <Card key={service.id} className="bg-slate-900 border-slate-800">
                            <CardContent className="p-4 flex items-center justify-between">
                                <div className="flex items-center gap-3">
                                    {status === 'checking' ? (
                                        <Loader2 className="w-5 h-5 text-slate-400 animate-spin" />
                                    ) : status === 'operational' ? (
                                        <div className="w-3 h-3 rounded-full bg-green-500 shadow-[0_0_10px_rgba(34,197,94,0.5)]" />
                                    ) : (
                                        <div className="w-3 h-3 rounded-full bg-red-500 shadow-[0_0_10px_rgba(239,68,68,0.5)]" />
                                    )}
                                    <span className="font-medium text-white">{service.name}</span>
                                </div>
                                <Badge variant="outline" className={`capitalize ${status === 'operational' ? 'border-green-500/30 text-green-400' :
                                    status === 'checking' ? 'border-slate-600 text-slate-400' :
                                        'border-red-500/30 text-red-400'
                                    }`}>
                                    {status}
                                </Badge>
                            </CardContent>
                        </Card>
                    )
                })}
            </motion.div>

            {/* Placeholder for Real Metrics */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.2 }}
            >
                <Card className="bg-slate-900 border-slate-800 border-dashed">
                    <CardContent className="p-12 flex flex-col items-center justify-center text-center">
                        <Activity className="w-12 h-12 text-slate-600 mb-4" />
                        <h3 className="text-lg font-medium text-slate-300">Real-time Metrics Coming Soon</h3>
                        <p className="text-slate-500 mt-2 max-w-md">
                            We are currently integrating real-time analytics for user growth, simulation stats, and system performance.
                        </p>
                    </CardContent>
                </Card>
            </motion.div>
        </div>
    )
}

