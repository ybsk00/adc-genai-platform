/**
 * Design Session WebSocket Hook
 * 실시간 세션 상태 관리 및 WebSocket 연결
 */
import { useState, useEffect, useCallback, useRef } from 'react'

export interface AgentStatus {
  agent: string
  status: 'started' | 'completed' | 'error' | 'warning'
  message: string
  step: number
  timestamp: string
}

export interface ProgressUpdate {
  current_step: number
  total_steps: number
  current_agent: string
  progress_percent: number
}

export interface CandidateUpdate {
  candidates: Candidate[]
  is_partial: boolean
  count: number
}

export interface Candidate {
  rank: number
  smiles: string
  score: number
  metrics?: Record<string, number>
  is_masked?: boolean
}

export interface SessionCompleteEvent {
  status: 'completed' | 'failed' | 'manual_review'
  final_report?: Record<string, unknown>
}

// New event types for UI enhancements
export interface SharedStateSync {
  current_smiles: string
  calculated_metrics: Record<string, number>
  scaffolds?: ScaffoldInfo[]
}

export interface ScaffoldInfo {
  id: string
  type: string
  smiles: string
  linkedKnowledgeIds?: string[]
}

export interface HealerAction {
  action: 'detecting' | 'analyzing' | 'healing' | 'healed' | 'failed'
  errorType?: string
  errorMessage?: string
  fixExplanation?: string
  attempt?: number
  maxAttempts?: number
  originalCode?: string
  fixedCode?: string
  timestamp: string
}

export interface LibrarianReference {
  id: string
  pmid?: string
  title: string
  authors?: string
  relevanceScore: number
  summary?: string
  highlightedSentences?: Array<{
    text: string
    relevance: 'high' | 'medium' | 'low'
    scaffoldId?: string
  }>
  linkedScaffoldIds?: string[]
}

export interface GoldenSetRef {
  id: string
  drugName: string
  target: string
  indication?: string
  clinicalStatus: string
  similarity: number
}

export interface AuditorFeedback {
  isInRedesignLoop: boolean
  checkItems: Array<{
    id: string
    name: string
    category: 'lipinski' | 'pains' | 'stability' | 'toxicity' | 'dar' | 'other'
    status: 'pass' | 'fail' | 'warning' | 'pending'
    currentValue?: number | string
    targetValue?: number | string
    unit?: string
    reasoning?: string
  }>
  redesignRequest?: {
    requestedBy: 'auditor'
    reason: string
    failedChecks: string[]
    suggestions: string[]
    targetAgent: 'alchemist' | 'coder'
    iteration: number
    maxIterations: number
    timestamp: string
  }
}

export interface DigitalSealInfo {
  sessionId: string
  recordHash: string
  previousHash?: string
  chainHash?: string
  timestamp: string
  isVerified: boolean
}

interface UseDesignSessionOptions {
  sessionId: string
  onAgentStatus?: (status: AgentStatus) => void
  onProgress?: (progress: ProgressUpdate) => void
  onCandidates?: (candidates: CandidateUpdate) => void
  onComplete?: (event: SessionCompleteEvent) => void
  onError?: (error: Error) => void
  // New event handlers
  onSharedStateSync?: (state: SharedStateSync) => void
  onHealerAction?: (action: HealerAction) => void
  onLibrarianReferences?: (refs: LibrarianReference[], goldenSets: GoldenSetRef[]) => void
  onAuditorFeedback?: (feedback: AuditorFeedback) => void
  onDigitalSeal?: (seal: DigitalSealInfo) => void
}

interface UseDesignSessionReturn {
  isConnected: boolean
  agentLogs: AgentStatus[]
  currentAgent: string | null
  currentStep: number
  progress: number
  candidates: Candidate[]
  sessionStatus: string | null
  // New state
  currentSmiles: string | null
  calculatedMetrics: Record<string, number> | null
  scaffolds: ScaffoldInfo[]
  healerAction: HealerAction | null
  isHealerActive: boolean
  literatureRefs: LibrarianReference[]
  goldenSetRefs: GoldenSetRef[]
  auditFeedback: AuditorFeedback | null
  digitalSeal: DigitalSealInfo | null
  // Methods
  connect: () => void
  disconnect: () => void
}

const WS_BASE_URL = import.meta.env.VITE_WS_URL || 'ws://localhost:8000'

export function useDesignSession({
  sessionId,
  onAgentStatus,
  onProgress,
  onCandidates,
  onComplete,
  onError,
  onSharedStateSync,
  onHealerAction,
  onLibrarianReferences,
  onAuditorFeedback,
  onDigitalSeal
}: UseDesignSessionOptions): UseDesignSessionReturn {
  const [isConnected, setIsConnected] = useState(false)
  const [agentLogs, setAgentLogs] = useState<AgentStatus[]>([])
  const [currentAgent, setCurrentAgent] = useState<string | null>(null)
  const [currentStep, setCurrentStep] = useState(0)
  const [progress, setProgress] = useState(0)
  const [candidates, setCandidates] = useState<Candidate[]>([])
  const [sessionStatus, setSessionStatus] = useState<string | null>(null)

  // New state for enhanced features
  const [currentSmiles, setCurrentSmiles] = useState<string | null>(null)
  const [calculatedMetrics, setCalculatedMetrics] = useState<Record<string, number> | null>(null)
  const [scaffolds, setScaffolds] = useState<ScaffoldInfo[]>([])
  const [healerAction, setHealerAction] = useState<HealerAction | null>(null)
  const [isHealerActive, setIsHealerActive] = useState(false)
  const [literatureRefs, setLiteratureRefs] = useState<LibrarianReference[]>([])
  const [goldenSetRefs, setGoldenSetRefs] = useState<GoldenSetRef[]>([])
  const [auditFeedback, setAuditFeedback] = useState<AuditorFeedback | null>(null)
  const [digitalSeal, setDigitalSeal] = useState<DigitalSealInfo | null>(null)

  const wsRef = useRef<WebSocket | null>(null)
  const reconnectTimeoutRef = useRef<ReturnType<typeof setTimeout> | null>(null)
  const pingIntervalRef = useRef<ReturnType<typeof setInterval> | null>(null)

  const connect = useCallback(() => {
    if (wsRef.current?.readyState === WebSocket.OPEN) {
      return
    }

    const ws = new WebSocket(`${WS_BASE_URL}/api/design/ws/${sessionId}`)
    wsRef.current = ws

    ws.onopen = () => {
      console.log('[ws] Connected to design session:', sessionId)
      setIsConnected(true)

      // Ping interval to keep connection alive
      pingIntervalRef.current = setInterval(() => {
        if (ws.readyState === WebSocket.OPEN) {
          ws.send('ping')
        }
      }, 30000)
    }

    ws.onmessage = (event) => {
      try {
        const data = JSON.parse(event.data)

        switch (data.type) {
          case 'initial_state':
            setCurrentAgent(data.current_agent)
            setCurrentStep(data.current_step || 0)
            setSessionStatus(data.status)
            break

          case 'agent_status':
            const status: AgentStatus = {
              agent: data.agent,
              status: data.status,
              message: data.message,
              step: data.step,
              timestamp: data.timestamp
            }
            setAgentLogs(prev => [...prev, status])
            setCurrentAgent(data.agent)
            setCurrentStep(data.step)
            onAgentStatus?.(status)
            break

          case 'progress':
            setCurrentStep(data.current_step)
            setProgress(data.progress_percent)
            setCurrentAgent(data.current_agent)
            onProgress?.(data)
            break

          case 'candidates_update':
            setCandidates(data.candidates)
            onCandidates?.(data)
            break

          case 'session_complete':
            setSessionStatus(data.status)
            onComplete?.(data)
            break

          // New event handlers for UI enhancements
          case 'shared_state_sync':
            setCurrentSmiles(data.current_smiles)
            setCalculatedMetrics(data.calculated_metrics)
            if (data.scaffolds) {
              setScaffolds(data.scaffolds)
            }
            onSharedStateSync?.(data)
            break

          case 'healer_action':
            const healerData: HealerAction = {
              action: data.action,
              errorType: data.error_type,
              errorMessage: data.error_message,
              fixExplanation: data.fix_explanation,
              attempt: data.attempt,
              maxAttempts: data.max_attempts || 3,
              originalCode: data.original_code,
              fixedCode: data.fixed_code,
              timestamp: data.timestamp || new Date().toISOString()
            }
            setHealerAction(healerData)
            setIsHealerActive(data.action !== 'healed' && data.action !== 'failed')

            // Add to agent logs for visibility
            setAgentLogs(prev => [...prev, {
              agent: 'healer',
              status: data.action === 'healed' ? 'completed' :
                     data.action === 'failed' ? 'error' : 'warning',
              message: data.action === 'detecting' ? '⚠️ Error Detected' :
                      data.action === 'analyzing' ? '🔍 Analyzing Error...' :
                      data.action === 'healing' ? '🛠️ Healing Code...' :
                      data.action === 'healed' ? '✅ Code Healed Successfully' :
                      '❌ Healing Failed',
              step: data.attempt || 0,
              timestamp: data.timestamp || new Date().toISOString()
            }])
            onHealerAction?.(healerData)
            break

          case 'librarian_references':
            const refs: LibrarianReference[] = (data.references || []).map((r: Record<string, unknown>) => ({
              id: r.id as string,
              pmid: r.pmid as string | undefined,
              title: r.title as string,
              authors: r.authors as string | undefined,
              relevanceScore: r.relevance_score as number || r.score as number || 0,
              summary: r.summary as string | undefined,
              highlightedSentences: r.highlighted_sentences as LibrarianReference['highlightedSentences'],
              linkedScaffoldIds: r.linked_scaffold_ids as string[]
            }))
            const goldenSets: GoldenSetRef[] = (data.golden_set_refs || []).map((g: Record<string, unknown>) => ({
              id: g.id as string,
              drugName: g.drug_name as string || g.name as string,
              target: g.target as string,
              indication: g.indication as string | undefined,
              clinicalStatus: g.clinical_status as string || g.status as string,
              similarity: g.similarity as number || 0
            }))
            setLiteratureRefs(refs)
            setGoldenSetRefs(goldenSets)
            onLibrarianReferences?.(refs, goldenSets)
            break

          case 'auditor_feedback':
            const feedback: AuditorFeedback = {
              isInRedesignLoop: data.is_in_redesign_loop || false,
              checkItems: (data.check_items || []).map((item: Record<string, unknown>) => ({
                id: item.id as string,
                name: item.name as string,
                category: (item.category || 'other') as AuditorFeedback['checkItems'][0]['category'],
                status: item.status as AuditorFeedback['checkItems'][0]['status'],
                currentValue: item.current_value as number | string | undefined,
                targetValue: item.target_value as number | string | undefined,
                unit: item.unit as string | undefined,
                reasoning: item.reasoning as string | undefined
              })),
              redesignRequest: data.redesign_request ? {
                requestedBy: 'auditor' as const,
                reason: data.redesign_request.reason,
                failedChecks: data.redesign_request.failed_checks || [],
                suggestions: data.redesign_request.suggestions || [],
                targetAgent: data.redesign_request.target_agent,
                iteration: data.redesign_request.iteration,
                maxIterations: data.redesign_request.max_iterations || 3,
                timestamp: data.redesign_request.timestamp || new Date().toISOString()
              } : undefined
            }
            setAuditFeedback(feedback)
            onAuditorFeedback?.(feedback)
            break

          case 'digital_seal':
            const seal: DigitalSealInfo = {
              sessionId: data.session_id,
              recordHash: data.record_hash,
              previousHash: data.previous_hash,
              chainHash: data.chain_hash,
              timestamp: data.timestamp,
              isVerified: data.is_verified !== false
            }
            setDigitalSeal(seal)
            onDigitalSeal?.(seal)
            break
        }
      } catch (e) {
        // Handle pong or non-JSON messages
        if (event.data !== 'pong') {
          console.warn('[ws] Failed to parse message:', e)
        }
      }
    }

    ws.onerror = (error) => {
      console.error('[ws] WebSocket error:', error)
      onError?.(new Error('WebSocket connection error'))
    }

    ws.onclose = (event) => {
      console.log('[ws] Connection closed:', event.code, event.reason)
      setIsConnected(false)

      // Clear ping interval
      if (pingIntervalRef.current) {
        clearInterval(pingIntervalRef.current)
      }

      // Auto-reconnect on abnormal closure (not manual disconnect)
      if (event.code !== 1000 && event.code !== 1001) {
        reconnectTimeoutRef.current = setTimeout(() => {
          console.log('[ws] Attempting to reconnect...')
          connect()
        }, 3000)
      }
    }
  }, [sessionId, onAgentStatus, onProgress, onCandidates, onComplete, onError])

  const disconnect = useCallback(() => {
    if (reconnectTimeoutRef.current) {
      clearTimeout(reconnectTimeoutRef.current)
    }
    if (pingIntervalRef.current) {
      clearInterval(pingIntervalRef.current)
    }
    if (wsRef.current) {
      wsRef.current.close(1000, 'Manual disconnect')
      wsRef.current = null
    }
  }, [])

  // Auto-connect on mount
  useEffect(() => {
    if (sessionId) {
      connect()
    }

    return () => {
      disconnect()
    }
  }, [sessionId, connect, disconnect])

  return {
    isConnected,
    agentLogs,
    currentAgent,
    currentStep,
    progress,
    candidates,
    sessionStatus,
    // New state
    currentSmiles,
    calculatedMetrics,
    scaffolds,
    healerAction,
    isHealerActive,
    literatureRefs,
    goldenSetRefs,
    auditFeedback,
    digitalSeal,
    // Methods
    connect,
    disconnect
  }
}
