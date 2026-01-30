/**
 * Agent Log Panel
 * 실시간 에이전트 연산 로그 표시
 * 
 * FIXED: 세부 추론 과정 실시간 표시
 */

import { useEffect, useRef, useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import {
  Bot,
  BookOpen,
  FlaskConical,
  Code,
  Wrench,
  Shield,
  CheckCircle2,
  AlertCircle,
  Loader2,
  FileText,
  ExternalLink
} from 'lucide-react';
import { ScrollArea } from '@/components/ui/scroll-area';
import { Badge } from '@/components/ui/badge';
import { cn } from '@/lib/utils';

export interface AgentLog {
  timestamp: string;
  agent_name: string;
  step: number;
  status: string;
  message: string;
  reasoning?: string;
  data_source?: string;
  pmid_refs?: string[];
  confidence_score?: number;
}

interface AgentLogPanelProps {
  logs: AgentLog[];
  isRunning: boolean;
  className?: string;
}

const AGENT_CONFIG: Record<string, { icon: React.ReactNode; color: string; label: string }> = {
  librarian: {
    icon: <BookOpen className="h-4 w-4" />,
    color: 'text-blue-400',
    label: 'Librarian'
  },
  alchemist: {
    icon: <FlaskConical className="h-4 w-4" />,
    color: 'text-purple-400',
    label: 'Alchemist'
  },
  coder: {
    icon: <Code className="h-4 w-4" />,
    color: 'text-cyan-400',
    label: 'Coder'
  },
  healer: {
    icon: <Wrench className="h-4 w-4" />,
    color: 'text-amber-400',
    label: 'Healer'
  },
  auditor: {
    icon: <Shield className="h-4 w-4" />,
    color: 'text-emerald-400',
    label: 'Auditor'
  },
  virtual_trial: {
    icon: <FileText className="h-4 w-4" />,
    color: 'text-pink-400',
    label: 'Clinical'
  },
  orchestrator: {
    icon: <Bot className="h-4 w-4" />,
    color: 'text-slate-400',
    label: 'System'
  }
};

function getStatusIcon(status: string) {
  switch (status) {
    case 'completed':
      return <CheckCircle2 className="h-4 w-4 text-emerald-400" />;
    case 'error':
    case 'failed':
      return <AlertCircle className="h-4 w-4 text-red-400" />;
    case 'healing':
      return <Loader2 className="h-4 w-4 text-amber-400 animate-spin" />;
    case 'running':
    case 'started':
      return <Loader2 className="h-4 w-4 text-blue-400 animate-spin" />;
    default:
      return <div className="h-4 w-4 rounded-full bg-slate-600" />;
  }
}

function formatTimestamp(timestamp: string): string {
  try {
    const date = new Date(timestamp);
    return date.toLocaleTimeString('ko-KR', { 
      hour12: false,
      hour: '2-digit',
      minute: '2-digit',
      second: '2-digit'
    });
  } catch {
    return timestamp;
  }
}

export function AgentLogPanel({ logs, isRunning, className }: AgentLogPanelProps) {
  const scrollRef = useRef<HTMLDivElement>(null);
  const [expandedLog, setExpandedLog] = useState<number | null>(null);

  // Auto-scroll to bottom
  useEffect(() => {
    if (scrollRef.current) {
      scrollRef.current.scrollTop = scrollRef.current.scrollHeight;
    }
  }, [logs]);

  return (
    <div className={cn("bg-slate-900/50 border border-slate-700 rounded-lg overflow-hidden", className)}>
      {/* Header */}
      <div className="flex items-center justify-between px-4 py-3 border-b border-slate-700/50 bg-slate-800/30">
        <div className="flex items-center gap-2">
          <Bot className="h-5 w-5 text-violet-400" />
          <h3 className="text-white font-medium">Agent Execution Log</h3>
          {isRunning && (
            <Badge variant="outline" className="text-emerald-400 border-emerald-400/50 animate-pulse text-xs">
              LIVE
            </Badge>
          )}
        </div>
        <div className="text-xs text-slate-500">
          {logs.length} events
        </div>
      </div>

      {/* Log List */}
      <ScrollArea className="h-80" ref={scrollRef}>
        <div className="p-4 space-y-3">
          <AnimatePresence mode="popLayout">
            {logs.map((log, index) => {
              const config = AGENT_CONFIG[log.agent_name] || AGENT_CONFIG.orchestrator;
              const isExpanded = expandedLog === index;

              return (
                <motion.div
                  key={`${log.timestamp}-${index}`}
                  initial={{ opacity: 0, x: -20 }}
                  animate={{ opacity: 1, x: 0 }}
                  exit={{ opacity: 0, x: 20 }}
                  transition={{ duration: 0.2 }}
                  className={cn(
                    "group relative p-3 rounded-lg border transition-all cursor-pointer",
                    isExpanded 
                      ? "bg-slate-800/80 border-slate-600" 
                      : "bg-slate-800/40 border-slate-700/50 hover:border-slate-600"
                  )}
                  onClick={() => setExpandedLog(isExpanded ? null : index)}
                >
                  {/* Main Row */}
                  <div className="flex items-start gap-3">
                    {/* Status Icon */}
                    <div className="mt-0.5">
                      {getStatusIcon(log.status)}
                    </div>

                    {/* Agent Icon */}
                    <div className={cn("p-1.5 rounded-md bg-slate-700/50", config.color)}>
                      {config.icon}
                    </div>

                    {/* Content */}
                    <div className="flex-1 min-w-0">
                      <div className="flex items-center gap-2 flex-wrap">
                        <span className={cn("font-medium text-sm", config.color)}>
                          {config.label}
                        </span>
                        <Badge variant="secondary" className="text-xs bg-slate-700 text-slate-300">
                          Step {log.step}
                        </Badge>
                        {log.confidence_score !== undefined && (
                          <Badge 
                            variant="outline" 
                            className={cn(
                              "text-xs",
                              log.confidence_score > 0.8 ? 'text-emerald-400' :
                              log.confidence_score > 0.5 ? 'text-amber-400' : 'text-red-400'
                            )}
                          >
                            {(log.confidence_score * 100).toFixed(0)}%
                          </Badge>
                        )}
                        <span className="text-xs text-slate-500 ml-auto">
                          {formatTimestamp(log.timestamp)}
                        </span>
                      </div>

                      <p className="text-slate-300 text-sm mt-1 leading-relaxed">
                        {log.message}
                      </p>

                      {/* Data Source */}
                      {log.data_source && (
                        <div className="flex items-center gap-1 mt-2">
                          <span className="text-xs text-slate-500">Source:</span>
                          <Badge variant="outline" className="text-xs text-slate-400 border-slate-600">
                            {log.data_source}
                          </Badge>
                        </div>
                      )}
                    </div>

                    {/* Expand Icon */}
                    {log.reasoning && (
                      <ExternalLink className={cn(
                        "h-4 w-4 text-slate-500 transition-transform",
                        isExpanded && "rotate-90"
                      )} />
                    )}
                  </div>

                  {/* Expanded Details */}
                  <AnimatePresence>
                    {isExpanded && log.reasoning && (
                      <motion.div
                        initial={{ height: 0, opacity: 0 }}
                        animate={{ height: 'auto', opacity: 1 }}
                        exit={{ height: 0, opacity: 0 }}
                        transition={{ duration: 0.2 }}
                        className="overflow-hidden"
                      >
                        <div className="mt-3 pt-3 border-t border-slate-700/50 space-y-2">
                          {/* Reasoning */}
                          <div>
                            <span className="text-xs text-slate-500 uppercase tracking-wider">Reasoning</span>
                            <p className="text-sm text-slate-300 mt-1 bg-slate-900/50 p-2 rounded">
                              {log.reasoning}
                            </p>
                          </div>

                          {/* PMIDs */}
                          {log.pmid_refs && log.pmid_refs.length > 0 && (
                            <div>
                              <span className="text-xs text-slate-500 uppercase tracking-wider">
                                References (PMID)
                              </span>
                              <div className="flex flex-wrap gap-1 mt-1">
                                {log.pmid_refs.map((pmid, i) => (
                                  <Badge 
                                    key={i} 
                                    variant="outline" 
                                    className="text-xs text-blue-400 border-blue-400/30 cursor-pointer hover:bg-blue-400/10"
                                    onClick={(e) => {
                                      e.stopPropagation();
                                      window.open(`https://pubmed.ncbi.nlm.nih.gov/${pmid}/`, '_blank');
                                    }}
                                  >
                                    {pmid}
                                  </Badge>
                                ))}
                              </div>
                            </div>
                          )}
                        </div>
                      </motion.div>
                    )}
                  </AnimatePresence>
                </motion.div>
              );
            })}
          </AnimatePresence>

          {logs.length === 0 && (
            <div className="text-center py-8 text-slate-500">
              <Bot className="h-12 w-12 mx-auto mb-3 opacity-30" />
              <p>Waiting for agent execution...</p>
            </div>
          )}
        </div>
      </ScrollArea>
    </div>
  );
}

export default AgentLogPanel;
