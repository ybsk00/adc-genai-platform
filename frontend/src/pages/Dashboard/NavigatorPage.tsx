/**
 * One-Click ADC Navigator Page
 * 질환명 입력만으로 최적의 ADC 설계를 자동으로 수행
 */

import { useState, useCallback, useRef, useEffect } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import { Rocket, History, Trash2, RefreshCw } from 'lucide-react';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { toast } from 'sonner';
import {
  SmartDiseaseInput,
  NavigatorPipeline,
  NavigatorResult,
  type AgentLog,
} from '@/components/navigator';
import { AntibodySelectionModal } from '@/components/navigator/AntibodySelectionModal';
import { API_BASE_URL } from '@/lib/api';

// ============================================================================
// Types
// ============================================================================

interface StepStatus {
  status: 'pending' | 'running' | 'completed' | 'error';
  message?: string;
}

interface AntibodyCandidate {
  id: string;
  name: string;
  target_protein: string;
  clinical_score: number;
  match_confidence: number;
}

interface GoldenCombination {
  antibody: string;
  linker: {
    type: string;
    smiles: string;
  };
  payload: {
    class: string;
    smiles: string;
  };
  dar: number;
  historical_orr: number | null;
}

interface VirtualTrial {
  predicted_orr: number;
  predicted_pfs_months: number;
  predicted_os_months: number;
  pk_data: Array<{ time_hours: number; concentration: number; free_payload: number }>;
  tumor_data: Array<{ day: number; treated: number; control: number }>;
  confidence: number;
}

interface NavigatorSessionResult {
  session_id: string;
  disease_name: string;
  target_protein: string | null;
  antibody_candidates: AntibodyCandidate[];
  golden_combination: GoldenCombination | null;
  physics_verified: boolean;
  virtual_trial: VirtualTrial | null;
  execution_time_seconds: number;
  warnings?: string[];  // FIXED: Add warnings field
  data_quality_score?: number;
}

interface RecentSession {
  id: string;
  disease_name: string;
  status: string;
  predicted_orr: number | null;
  created_at: string;
}

// ============================================================================
// Component
// ============================================================================

export function NavigatorPage() {
  // State
  const [isRunning, setIsRunning] = useState(false);
  const [currentStep, setCurrentStep] = useState(0);
  const [diseaseName, setDiseaseName] = useState('');
  const [steps, setSteps] = useState<StepStatus[]>([
    { status: 'pending' },
    { status: 'pending' },
    { status: 'pending' },
    { status: 'pending' },
    { status: 'pending' },
  ]);
  const [result, setResult] = useState<NavigatorSessionResult | null>(null);
  const [recentSessions, setRecentSessions] = useState<RecentSession[]>([]);
  const [showHistory, setShowHistory] = useState(false);
  const [agentLogs, setAgentLogs] = useState<AgentLog[]>([]);  // FIXED: 실시간 에이전트 로그
  
  // Selection Gate State
  const [showSelectionModal, setShowSelectionModal] = useState(false);
  const [selectionCandidates, setSelectionCandidates] = useState<AntibodyCandidate[]>([]);
  const [currentSessionId, setCurrentSessionId] = useState<string | null>(null);

  const eventSourceRef = useRef<EventSource | null>(null);

  // Load recent sessions
  useEffect(() => {
    loadRecentSessions();
  }, []);

  const loadRecentSessions = async () => {
    try {
      const response = await fetch(`${API_BASE_URL}/api/design/navigator/sessions?limit=5`);
      if (response.ok) {
        const data = await response.json();
        setRecentSessions(data.sessions || []);
      }
    } catch (error) {
      console.error('Failed to load recent sessions:', error);
    }
  };

  // Reset state
  const resetState = useCallback(() => {
    setCurrentStep(0);
    setSteps([
      { status: 'pending' },
      { status: 'pending' },
      { status: 'pending' },
      { status: 'pending' },
      { status: 'pending' },
    ]);
    setResult(null);
    setAgentLogs([]);  // FIXED
    setShowSelectionModal(false);
    setSelectionCandidates([]);
    setCurrentSessionId(null);
  }, []);

  // Handle Antibody Selection
  const handleAntibodySelect = useCallback(async (antibodyId: string) => {
    if (!currentSessionId) return;
    
    setShowSelectionModal(false);
    // Keep isRunning true as we resume
    
    try {
        // Resume Navigator via POST with selected_antibody_id
        await fetch(`${API_BASE_URL}/api/design/navigator/run`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ 
                disease_name: diseaseName,
                session_id: currentSessionId,
                selected_antibody_id: antibodyId
            }),
        });
        
        // SSE will pick up the progress as we didn't close it (or we can reconnect)
        // Ideally, the backend updates the status to 'running' and the existing SSE loop continues.
        
    } catch (error) {
        console.error('Failed to select antibody:', error);
        toast.error('Failed to resume design');
    }
  }, [currentSessionId, diseaseName]);

  // Handle Navigator run with WebSocket streaming
  const handleNavigate = useCallback(async (disease: string, targetProtein?: string) => {
    resetState();
    setDiseaseName(disease);
    setIsRunning(true);

    // Close existing EventSource
    if (eventSourceRef.current) {
      eventSourceRef.current.close();
    }

    try {
      // Start Navigator via POST (disease + target)
      const response = await fetch(`${API_BASE_URL}/api/design/navigator/run`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          disease_name: disease,
          target_protein: targetProtein || null,
        }),
      });

      if (!response.ok) {
        throw new Error('Failed to start Navigator');
      }

      const { session_id } = await response.json();
      setCurrentSessionId(session_id);

      // Connect to SSE stream
      const eventSource = new EventSource(
        `${API_BASE_URL}/api/design/navigator/stream/${session_id}`
      );
      eventSourceRef.current = eventSource;

      eventSource.onmessage = (event) => {
        try {
          const data = JSON.parse(event.data);

          switch (data.type) {
            // ... (existing cases) ...
            
            // [NEW] Handle Waiting for Selection (Backend might send this or we poll)
            // Ideally backend sends a specific event, but if it relies on polling state change:
            
            case 'step_start':
              setCurrentStep(data.step);
              setSteps((prev) => {
                const newSteps = [...prev];
                if (newSteps[data.step - 1]) {
                  newSteps[data.step - 1] = {
                    status: 'running',
                    message: data.message,
                  };
                }
                return newSteps;
              });
              break;

            case 'step_complete':
              setSteps((prev) => {
                const newSteps = [...prev];
                if (newSteps[data.step - 1]) {
                  newSteps[data.step - 1] = {
                    status: 'completed',
                    message: data.message,
                  };
                }
                return newSteps;
              });
              break;

            case 'step_error':
              setSteps((prev) => {
                const newSteps = [...prev];
                if (newSteps[data.step - 1]) {
                  newSteps[data.step - 1] = {
                    status: 'error',
                    message: data.message,
                  };
                }
                return newSteps;
              });
              toast.error(`Step ${data.step} failed: ${data.message}`);
              break;

            // FIXED: 실시간 에이전트 로그 수신
            case 'agent_log':
              if (data.agent && data.message) {
                setAgentLogs((prev) => [...prev, {
                  timestamp: data.timestamp || new Date().toISOString(),
                  agent_name: data.agent,
                  step: data.step || currentStep,
                  status: data.status || 'running',
                  message: data.message,
                  reasoning: data.reasoning,
                  data_source: data.data_source,
                  pmid_refs: data.pmid_refs,
                  confidence_score: data.confidence
                }]);
              }
              break;

            case 'complete':
              setResult(data.result);
              // FIXED: 완료 시 로그 동기화
              if (data.result?.agent_logs) {
                setAgentLogs(data.result.agent_logs);
              }
              setIsRunning(false);
              eventSource.close();
              toast.success('ADC design completed!');
              loadRecentSessions();
              break;

            case 'error':
              toast.error(data.message || 'Navigator failed');
              setIsRunning(false);
              eventSource.close();
              break;
          }
        } catch (e) {
          console.error('Failed to parse SSE message:', e);
        }
      };

      eventSource.onerror = () => {
        // Don't close immediately, might be temporary
        // But if persistent, fall back to poll
        // eventSource.close();
        pollForResult(session_id);
      };

    } catch (error) {
      console.error('Navigator error:', error);
      toast.error('Failed to start Navigator');
      setIsRunning(false);
    }
  }, [resetState]);

  // Fallback polling if SSE fails (and for handling waiting state if SSE doesn't send it)
  const pollForResult = useCallback(async (sessionId: string) => {
    const maxAttempts = 120; // Increased
    let attempts = 0;

    const poll = async () => {
      try {
        const response = await fetch(`${API_BASE_URL}/api/design/navigator/session/${sessionId}`);
        if (response.ok) {
          const data = await response.json();

          // [NEW] Handle Waiting for Selection
          if (data.status === 'waiting_for_selection') {
             if (!showSelectionModal) {
                 // Fetch candidates to show
                 // Ideally they are in data, otherwise we rely on what we have or fetch again
                 // NavigatorSession table has 'antibody_candidates'
                 if (data.antibody_candidates && data.antibody_candidates.length > 0) {
                     setSelectionCandidates(data.antibody_candidates);
                     setShowSelectionModal(true);
                     // Pause polling or slow it down? 
                     // We can continue polling to detect when status changes back to running
                 }
             }
             // Update step to completed for step 1
             setSteps((prev) => {
                const newSteps = [...prev];
                newSteps[0] = { status: 'completed', message: 'Candidates found' };
                return newSteps;
             });
          }
          
          if (data.status === 'running') {
              // Close modal if open (means user selected)
              if (showSelectionModal) setShowSelectionModal(false);
          }

          if (data.status === 'completed') {
            setResult(data);
            setIsRunning(false);
            setSteps([
              { status: 'completed' },
              { status: 'completed' },
              { status: 'completed' },
              { status: 'completed' },
              { status: 'completed' },
            ]);
            toast.success('ADC design completed!');
            loadRecentSessions();
            return;
          }

          if (data.status === 'failed') {
            toast.error(data.error_message || 'Navigator failed');
            setIsRunning(false);
            return;
          }

          // Update step progress
          if (data.current_step) {
            setCurrentStep(data.current_step);
            setSteps((prev) => {
              const newSteps = [...prev];
              for (let i = 0; i < data.current_step; i++) {
                if (newSteps[i].status !== 'completed') {
                    newSteps[i] = { status: 'completed' };
                }
              }
              if (newSteps[data.current_step - 1]) {
                newSteps[data.current_step - 1] = { status: 'running' };
              }
              // Special case: if waiting, step 1 is done, step 2 not started
              if (data.status === 'waiting_for_selection') {
                  newSteps[0] = { status: 'completed' };
                  newSteps[1] = { status: 'pending' }; // Waiting
              }
              return newSteps;
            });
          }
        }
      } catch (e) {
        console.error('Poll error:', e);
      }

      attempts++;
      // Stop polling if completed or failed or stopped manually
      // Continue polling if waiting_for_selection to catch resume
      if (attempts < maxAttempts && (isRunning || showSelectionModal)) {
        setTimeout(poll, 2000);
      } else if (attempts >= maxAttempts) {
        if (isRunning) {
            toast.error('Navigator timeout');
            setIsRunning(false);
        }
      }
    };

    poll();
  }, [isRunning, showSelectionModal]);

  // Load a previous session
  const handleLoadSession = useCallback(async (sessionId: string) => {
    try {
      const response = await fetch(`${API_BASE_URL}/api/design/navigator/session/${sessionId}`);
      if (response.ok) {
        const data = await response.json();
        setDiseaseName(data.disease_name);
        setResult(data);
        setSteps([
          { status: 'completed' },
          { status: 'completed' },
          { status: 'completed' },
          { status: 'completed' },
          { status: 'completed' },
        ]);
        setShowHistory(false);
      }
    } catch (error) {
      toast.error('Failed to load session');
    }
  }, []);

  // Handle new design
  const handleNewDesign = useCallback(() => {
    resetState();
    setDiseaseName('');
  }, [resetState]);

  return (
    <div className="space-y-8">
      {/* Page Header */}
      <div className="text-center">
        <motion.div
          initial={{ opacity: 0, y: -20 }}
          animate={{ opacity: 1, y: 0 }}
          className="inline-flex items-center gap-3 mb-4"
        >
          <div className="p-3 bg-gradient-to-r from-violet-600 to-indigo-600 rounded-xl">
            <Rocket className="h-8 w-8 text-white" />
          </div>
          <div className="text-left">
            <h1 className="text-3xl font-bold text-white">One-Click ADC Navigator</h1>
            <p className="text-slate-400">
              Enter a disease name and let AI design the optimal ADC
            </p>
          </div>
        </motion.div>
      </div>

      {/* Selection Modal */}
      <AntibodySelectionModal
        open={showSelectionModal}
        candidates={selectionCandidates}
        onSelect={handleAntibodySelect}
        isLoading={isRunning && !showSelectionModal} // Loading if running but modal closed? No, simple logic
      />

      {/* Main Content */}
      <AnimatePresence mode="wait">
        {!isRunning && !result ? (
          // Initial State - Disease Input
          <motion.div
            key="input"
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
            className="space-y-8"
          >
            {/* Disease Input */}
            <div className="py-8">
              <SmartDiseaseInput
                onSubmit={handleNavigate}
                isLoading={isRunning}
                disabled={isRunning}
              />
            </div>

            {/* Recent Sessions */}
            {recentSessions.length > 0 && (
              <Card className="bg-slate-900/50 border-slate-700">
                <CardHeader className="flex flex-row items-center justify-between">
                  <div>
                    <CardTitle className="text-white flex items-center gap-2">
                      <History className="h-5 w-5 text-slate-400" />
                      Recent Designs
                    </CardTitle>
                    <CardDescription>Your previous Navigator sessions</CardDescription>
                  </div>
                  <Button
                    variant="ghost"
                    size="sm"
                    onClick={() => setShowHistory(!showHistory)}
                    className="text-slate-400 hover:text-white"
                  >
                    {showHistory ? 'Hide' : 'Show All'}
                  </Button>
                </CardHeader>
                <CardContent>
                  <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-3">
                    {(showHistory ? recentSessions : recentSessions.slice(0, 3)).map(
                      (session) => (
                        <motion.div
                          key={session.id}
                          whileHover={{ scale: 1.02 }}
                          className="p-4 bg-slate-800/50 rounded-lg border border-slate-700 cursor-pointer hover:border-violet-500/50 transition-colors"
                          onClick={() => handleLoadSession(session.id)}
                        >
                          <div className="flex items-start justify-between">
                            <div>
                              <p className="font-medium text-white">
                                {session.disease_name}
                              </p>
                              <p className="text-xs text-slate-500 mt-1">
                                {new Date(session.created_at).toLocaleDateString()}
                              </p>
                            </div>
                            <Badge
                              variant={
                                session.status === 'completed' ? 'default' : 'secondary'
                              }
                              className={
                                session.status === 'completed'
                                  ? 'bg-emerald-500/20 text-emerald-400'
                                  : 'bg-slate-700 text-slate-400'
                              }
                            >
                              {session.status === 'completed' && session.predicted_orr != null
                                ? `ORR ${Number(session.predicted_orr).toFixed(1)}%`
                                : session.status}
                            </Badge>
                          </div>
                        </motion.div>
                      )
                    )}
                  </div>
                </CardContent>
              </Card>
            )}

            {/* Feature Highlights */}
            <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
              {[
                {
                  title: '6 AI Agents',
                  description: 'Librarian, Alchemist, Coder, Auditor, Clinical, and PM agents collaborate',
                  gradient: 'from-violet-500 to-purple-500',
                },
                {
                  title: '7,600+ Compounds',
                  description: 'Analyzing Golden Set data for optimal combinations',
                  gradient: 'from-blue-500 to-cyan-500',
                },
                {
                  title: 'Virtual Trial',
                  description: 'PK/PD simulation with predicted ORR, PFS, OS',
                  gradient: 'from-emerald-500 to-teal-500',
                },
              ].map((feature, index) => (
                <motion.div
                  key={feature.title}
                  initial={{ opacity: 0, y: 20 }}
                  animate={{ opacity: 1, y: 0 }}
                  transition={{ delay: index * 0.1 }}
                  className="p-4 bg-slate-900/50 border border-slate-700 rounded-lg"
                >
                  <div
                    className={`w-10 h-10 rounded-lg bg-gradient-to-r ${feature.gradient} flex items-center justify-center mb-3`}
                  >
                    <span className="text-white font-bold">{index + 1}</span>
                  </div>
                  <h3 className="text-white font-semibold">{feature.title}</h3>
                  <p className="text-sm text-slate-400 mt-1">{feature.description}</p>
                </motion.div>
              ))}
            </div>
          </motion.div>
        ) : isRunning ? (
          // Running State - Pipeline View
          <motion.div
            key="running"
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
          >
            <NavigatorPipeline
              currentStep={currentStep}
              steps={steps}
              diseaseName={diseaseName}
              agentLogs={agentLogs}  // FIXED: 실시간 로그 전달
            />
          </motion.div>
        ) : result ? (
          // Complete State - Results View
          <motion.div
            key="result"
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            exit={{ opacity: 0 }}
            className="space-y-6"
          >
            {/* New Design Button */}
            <div className="flex justify-end">
              <Button
                variant="outline"
                onClick={handleNewDesign}
                className="border-slate-700 hover:bg-slate-800"
              >
                <RefreshCw className="mr-2 h-4 w-4" />
                New Design
              </Button>
            </div>

            <NavigatorResult
              sessionId={result.session_id}
              diseaseName={result.disease_name}
              targetProtein={result.target_protein}
              antibodyCandidates={result.antibody_candidates}
              goldenCombination={result.golden_combination}
              physicsVerified={result.physics_verified}
              virtualTrial={result.virtual_trial}
              executionTime={result.execution_time_seconds}
              warnings={result.warnings}  // FIXED
              dataQualityScore={result.data_quality_score}  // FIXED
              onExportReport={() => toast.info('Export feature coming soon')}
              onShare={() => toast.info('Share feature coming soon')}
            />
          </motion.div>
        ) : null}
      </AnimatePresence>
    </div>
  );
}

export default NavigatorPage;
