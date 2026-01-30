/**
 * Navigator Pipeline Component
 * One-Click ADC Navigator - Wizard-Style Pipeline Progress
 *
 * 4-Phase Wizard + 5-Step Pipeline 진행 상황 표시
 */

import { useState, useEffect } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import {
  CheckCircle2,
  Circle,
  Loader2,
  AlertCircle,
  BookOpen,
  FlaskConical,
  Calculator,
  Shield,
  Activity,
  ChevronRight,
  Search,
  Beaker,
  ClipboardCheck,
  BarChart3,
} from 'lucide-react';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Progress } from '@/components/ui/progress';
import { cn } from '@/lib/utils';

interface PipelineStep {
  id: number;
  name: string;
  description: string;
  agent: string;
  icon: React.ReactNode;
  phase: number;
  status: 'pending' | 'running' | 'completed' | 'error';
  message?: string;
}

import { AgentLogPanel, type AgentLog } from './AgentLogPanel';

interface NavigatorPipelineProps {
  currentStep: number;
  steps: Array<{
    status: 'pending' | 'running' | 'completed' | 'error';
    message?: string;
  }>;
  diseaseName: string;
  agentLogs?: AgentLog[];
  onComplete?: () => void;
}

// 4-Phase Wizard 정의
const WIZARD_PHASES = [
  { id: 1, name: 'Discovery', icon: <Search className="h-4 w-4" />, steps: [1] },
  { id: 2, name: 'Design', icon: <Beaker className="h-4 w-4" />, steps: [2, 3] },
  { id: 3, name: 'Validate', icon: <ClipboardCheck className="h-4 w-4" />, steps: [4] },
  { id: 4, name: 'Simulate', icon: <BarChart3 className="h-4 w-4" />, steps: [5] },
];

const DEFAULT_STEPS: Omit<PipelineStep, 'status' | 'message'>[] = [
  {
    id: 1,
    name: 'Target & Antibody Match',
    description: 'Finding optimal antibodies for the disease',
    agent: 'Librarian',
    icon: <BookOpen className="h-5 w-5" />,
    phase: 1,
  },
  {
    id: 2,
    name: 'Golden Combination',
    description: 'Selecting best linker-payload combination',
    agent: 'Alchemist',
    icon: <FlaskConical className="h-5 w-5" />,
    phase: 2,
  },
  {
    id: 3,
    name: 'Property Calculation',
    description: 'Computing molecular properties',
    agent: 'Coder + Healer',
    icon: <Calculator className="h-5 w-5" />,
    phase: 2,
  },
  {
    id: 4,
    name: 'Physical Validation',
    description: 'Verifying structural feasibility',
    agent: 'Auditor',
    icon: <Shield className="h-5 w-5" />,
    phase: 3,
  },
  {
    id: 5,
    name: 'Virtual Trial',
    description: 'Running clinical simulation',
    agent: 'Clinical',
    icon: <Activity className="h-5 w-5" />,
    phase: 4,
  },
];

const statusColors = {
  pending: 'text-slate-500 bg-slate-800',
  running: 'text-violet-400 bg-violet-900/50',
  completed: 'text-emerald-400 bg-emerald-900/50',
  error: 'text-red-400 bg-red-900/50',
};

function getPhaseStatus(
  phase: (typeof WIZARD_PHASES)[number],
  steps: Array<{ status: string }>
): 'pending' | 'running' | 'completed' | 'error' {
  const phaseStepStatuses = phase.steps.map((s) => steps[s - 1]?.status || 'pending');
  if (phaseStepStatuses.every((s) => s === 'completed')) return 'completed';
  if (phaseStepStatuses.some((s) => s === 'error')) return 'error';
  if (phaseStepStatuses.some((s) => s === 'running')) return 'running';
  return 'pending';
}

export function NavigatorPipeline({
  currentStep,
  steps,
  diseaseName,
  agentLogs = [],
  onComplete,
}: NavigatorPipelineProps) {
  const [progress, setProgress] = useState(0);

  useEffect(() => {
    const completedSteps = steps.filter((s) => s.status === 'completed').length;
    setProgress((completedSteps / DEFAULT_STEPS.length) * 100);

    if (completedSteps === DEFAULT_STEPS.length) {
      onComplete?.();
    }
  }, [steps, onComplete]);

  const pipelineSteps: PipelineStep[] = DEFAULT_STEPS.map((step, index) => ({
    ...step,
    status: steps[index]?.status || 'pending',
    message: steps[index]?.message,
  }));

  return (
    <Card className="bg-slate-900/50 border-slate-700 overflow-hidden">
      <CardHeader className="border-b border-slate-700/50 pb-6">
        <div className="flex items-center justify-between mb-4">
          <div>
            <CardTitle className="text-white flex items-center gap-2">
              ADC Design Wizard
              {currentStep > 0 && currentStep <= 5 && (
                <Badge
                  variant="outline"
                  className="text-violet-400 border-violet-400/50 animate-pulse"
                >
                  Running
                </Badge>
              )}
            </CardTitle>
            <p className="text-sm text-slate-500 mt-1">
              Designing ADC for:{' '}
              <span className="text-violet-400 font-medium">{diseaseName}</span>
            </p>
          </div>
          <div className="text-right">
            <p className="text-sm text-slate-500">Progress</p>
            <p className="text-2xl font-bold text-white">{Math.round(progress)}%</p>
          </div>
        </div>

        {/* Wizard Phase Indicator (horizontal stepper) */}
        <div className="flex items-center gap-0 mt-2">
          {WIZARD_PHASES.map((phase, idx) => {
            const phaseStatus = getPhaseStatus(phase, steps);
            const isActive = phaseStatus === 'running';
            const isCompleted = phaseStatus === 'completed';
            const isError = phaseStatus === 'error';

            return (
              <div key={phase.id} className="flex items-center flex-1">
                <div className="flex flex-col items-center flex-1">
                  {/* Phase circle */}
                  <motion.div
                    animate={isActive ? { scale: [1, 1.1, 1] } : {}}
                    transition={isActive ? { repeat: Infinity, duration: 1.5 } : {}}
                    className={cn(
                      'flex items-center justify-center w-10 h-10 rounded-full border-2 transition-all',
                      isCompleted && 'bg-emerald-500/20 border-emerald-500 text-emerald-400',
                      isActive && 'bg-violet-500/20 border-violet-500 text-violet-400',
                      isError && 'bg-red-500/20 border-red-500 text-red-400',
                      phaseStatus === 'pending' &&
                        'bg-slate-800 border-slate-600 text-slate-500'
                    )}
                  >
                    {isCompleted ? (
                      <CheckCircle2 className="h-5 w-5" />
                    ) : isActive ? (
                      <Loader2 className="h-5 w-5 animate-spin" />
                    ) : isError ? (
                      <AlertCircle className="h-5 w-5" />
                    ) : (
                      phase.icon
                    )}
                  </motion.div>
                  {/* Phase label */}
                  <span
                    className={cn(
                      'text-xs mt-1.5 font-medium',
                      isCompleted && 'text-emerald-400',
                      isActive && 'text-violet-400',
                      isError && 'text-red-400',
                      phaseStatus === 'pending' && 'text-slate-500'
                    )}
                  >
                    {phase.name}
                  </span>
                </div>
                {/* Connector line between phases */}
                {idx < WIZARD_PHASES.length - 1 && (
                  <div className="flex-shrink-0 w-8 h-0.5 -mt-5 mx-1">
                    <div
                      className={cn(
                        'h-full rounded-full transition-colors',
                        isCompleted ? 'bg-emerald-500/60' : 'bg-slate-700'
                      )}
                    />
                  </div>
                )}
              </div>
            );
          })}
        </div>

        <Progress value={progress} className="mt-5 h-2" />
      </CardHeader>

      {/* Step details */}
      <CardContent className="p-0">
        <div className="divide-y divide-slate-800">
          {pipelineSteps.map((step, index) => {
            const isActive = step.status === 'running';
            const isCompleted = step.status === 'completed';
            const isPending = step.status === 'pending';

            return (
              <motion.div
                key={step.id}
                initial={{ opacity: 0, x: -20 }}
                animate={{ opacity: 1, x: 0 }}
                transition={{ delay: index * 0.1 }}
                className={cn(
                  'flex items-center gap-4 p-4 transition-colors',
                  isActive && 'bg-violet-900/20',
                  isCompleted && 'bg-emerald-900/10'
                )}
              >
                {/* Step Number & Icon */}
                <div
                  className={cn(
                    'flex items-center justify-center w-10 h-10 rounded-full transition-colors',
                    statusColors[step.status]
                  )}
                >
                  {step.status === 'pending' ? (
                    <span className="text-lg font-bold">{step.id}</span>
                  ) : step.status === 'running' ? (
                    <Loader2 className="h-5 w-5 animate-spin" />
                  ) : step.status === 'completed' ? (
                    <CheckCircle2 className="h-5 w-5" />
                  ) : (
                    <AlertCircle className="h-5 w-5" />
                  )}
                </div>

                {/* Step Info */}
                <div className="flex-1">
                  <div className="flex items-center gap-2">
                    <h4
                      className={cn(
                        'font-medium transition-colors',
                        isPending && 'text-slate-500',
                        isActive && 'text-violet-300',
                        isCompleted && 'text-emerald-300',
                        step.status === 'error' && 'text-red-300'
                      )}
                    >
                      {step.name}
                    </h4>
                    <Badge
                      variant="secondary"
                      className={cn(
                        'text-xs',
                        isPending && 'bg-slate-800 text-slate-500',
                        isActive && 'bg-violet-900/50 text-violet-400',
                        isCompleted && 'bg-emerald-900/50 text-emerald-400'
                      )}
                    >
                      {step.agent}
                    </Badge>
                  </div>
                  <AnimatePresence mode="wait">
                    <motion.p
                      key={step.message || step.description}
                      initial={{ opacity: 0 }}
                      animate={{ opacity: 1 }}
                      exit={{ opacity: 0 }}
                      className={cn(
                        'text-sm mt-0.5',
                        isPending && 'text-slate-600',
                        isActive && 'text-slate-400',
                        isCompleted && 'text-slate-500'
                      )}
                    >
                      {step.message || step.description}
                    </motion.p>
                  </AnimatePresence>
                </div>

                {/* Agent Icon */}
                <div
                  className={cn(
                    'p-2 rounded-lg transition-colors',
                    isPending && 'bg-slate-800/50 text-slate-600',
                    isActive && 'bg-violet-900/30 text-violet-400',
                    isCompleted && 'bg-emerald-900/30 text-emerald-400'
                  )}
                >
                  {step.icon}
                </div>

                {/* Connector */}
                {index < pipelineSteps.length - 1 && (
                  <ChevronRight
                    className={cn(
                      'h-4 w-4',
                      isPending && 'text-slate-700',
                      isCompleted && 'text-emerald-600'
                    )}
                  />
                )}
              </motion.div>
            );
          })}
        </div>
      </CardContent>

      {/* Real-time agent log panel */}
      {agentLogs.length > 0 && (
        <div className="border-t border-slate-700/50">
          <AgentLogPanel
            logs={agentLogs}
            isRunning={currentStep > 0 && currentStep < 6}
            className="border-0 rounded-none"
          />
        </div>
      )}
    </Card>
  );
}

export default NavigatorPipeline;
