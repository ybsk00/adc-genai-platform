/**
 * Live Agent Console
 * 실시간 에이전트 실행 상태 표시 콘솔
 */
import { useRef, useEffect } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import { ScrollArea } from '@/components/ui/scroll-area'
import { motion, AnimatePresence } from 'framer-motion'
import {
  Terminal,
  CheckCircle2,
  XCircle,
  AlertTriangle,
  Loader2,
  Beaker,
  Code2,
  HeartPulse,
  Shield,
  BookOpen
} from 'lucide-react'
import type { AgentStatus } from '@/hooks/useDesignSession'

interface LiveAgentConsoleProps {
  logs: AgentStatus[]
  currentAgent: string | null
  isConnected: boolean
}

const AGENT_ICONS: Record<string, React.ReactNode> = {
  librarian: <BookOpen className="w-4 h-4" />,
  alchemist: <Beaker className="w-4 h-4" />,
  coder: <Code2 className="w-4 h-4" />,
  healer: <HeartPulse className="w-4 h-4" />,
  auditor: <Shield className="w-4 h-4" />,
  system: <Terminal className="w-4 h-4" />
}

const AGENT_COLORS: Record<string, string> = {
  librarian: 'text-amber-500 bg-amber-50 border-amber-200',
  alchemist: 'text-purple-500 bg-purple-50 border-purple-200',
  coder: 'text-blue-500 bg-blue-50 border-blue-200',
  healer: 'text-green-500 bg-green-50 border-green-200',
  auditor: 'text-red-500 bg-red-50 border-red-200',
  system: 'text-gray-500 bg-gray-50 border-gray-200'
}

const STATUS_ICONS: Record<string, React.ReactNode> = {
  started: <Loader2 className="w-3 h-3 animate-spin" />,
  completed: <CheckCircle2 className="w-3 h-3 text-green-500" />,
  error: <XCircle className="w-3 h-3 text-red-500" />,
  warning: <AlertTriangle className="w-3 h-3 text-yellow-500" />
}

export function LiveAgentConsole({ logs, currentAgent, isConnected }: LiveAgentConsoleProps) {
  const scrollRef = useRef<HTMLDivElement>(null)

  // Auto-scroll to bottom when new logs arrive
  useEffect(() => {
    if (scrollRef.current) {
      scrollRef.current.scrollTop = scrollRef.current.scrollHeight
    }
  }, [logs])

  const formatTime = (timestamp: string) => {
    const date = new Date(timestamp)
    return date.toLocaleTimeString('ko-KR', {
      hour: '2-digit',
      minute: '2-digit',
      second: '2-digit'
    })
  }

  return (
    <Card className="h-full flex flex-col">
      <CardHeader className="flex flex-row items-center justify-between pb-2 space-y-0">
        <CardTitle className="flex items-center gap-2 text-sm font-medium">
          <Terminal className="w-4 h-4" />
          Agent Console
        </CardTitle>
        <div className="flex items-center gap-2">
          {currentAgent && (
            <Badge variant="outline" className={`text-xs ${AGENT_COLORS[currentAgent] || ''}`}>
              {AGENT_ICONS[currentAgent]}
              <span className="ml-1 capitalize">{currentAgent}</span>
            </Badge>
          )}
          <span
            className={`w-2 h-2 rounded-full ${isConnected ? 'bg-green-500' : 'bg-red-500'}`}
            title={isConnected ? 'Connected' : 'Disconnected'}
          />
        </div>
      </CardHeader>
      <CardContent className="flex-1 overflow-hidden p-0">
        <ScrollArea className="h-[400px] p-4" ref={scrollRef}>
          <div className="space-y-2 font-mono text-xs">
            <AnimatePresence mode="popLayout">
              {logs.length === 0 ? (
                <motion.div
                  initial={{ opacity: 0 }}
                  animate={{ opacity: 1 }}
                  className="text-gray-400 text-center py-8"
                >
                  Waiting for agent activity...
                </motion.div>
              ) : (
                logs.map((log, index) => (
                  <motion.div
                    key={`${log.timestamp}-${index}`}
                    initial={{ opacity: 0, x: -10 }}
                    animate={{ opacity: 1, x: 0 }}
                    transition={{ duration: 0.2 }}
                    className={`flex items-start gap-2 p-2 rounded border ${
                      log.status === 'error'
                        ? 'bg-red-50 border-red-200'
                        : log.status === 'warning'
                        ? 'bg-yellow-50 border-yellow-200'
                        : 'bg-slate-50 border-slate-200'
                    }`}
                  >
                    <span className="text-gray-400 whitespace-nowrap">
                      [{formatTime(log.timestamp)}]
                    </span>
                    <span className={`flex items-center gap-1 whitespace-nowrap ${
                      AGENT_COLORS[log.agent]?.split(' ')[0] || 'text-gray-500'
                    }`}>
                      {AGENT_ICONS[log.agent] || AGENT_ICONS.system}
                      <span className="uppercase font-semibold">{log.agent}</span>
                    </span>
                    <span className="flex items-center gap-1">
                      {STATUS_ICONS[log.status]}
                    </span>
                    <span className="text-gray-700 break-words flex-1">
                      {log.message}
                    </span>
                  </motion.div>
                ))
              )}
            </AnimatePresence>
          </div>
        </ScrollArea>
      </CardContent>
    </Card>
  )
}
