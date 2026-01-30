/**
 * Navigator Pipeline Component
 * One-Click ADC Navigator - Live Pipeline Progress
 *
 * 6인 에이전트 협업 파이프라인 진행 상황 표시
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
  status: 'pending' | 'running' | 'completed' | 'error';
  message?: string;
}

import { AgentLogPanel, AgentLog } from './AgentLogPanel';

interface NavigatorPipelineProps {
  currentStep: number;
  steps: Array<{
    status: 'pending' | 'running' | 'completed' | 'error';
    message?: string;
  }>;
  diseaseName: string;
  agentLogs?: AgentLog[];  // FIXED: 실시간 에이전트 로그
  onComplete?: () => void;
}

const DEFAULT_STEPS: Omit<PipelineStep, 'status' | 'message'>[] = [
  {
    id: 1,
    name: 'Target & Antibody Match',
    description: 'Finding optimal antibodies for the disease',
    agent: 'Librarian',
    icon: <BookOpen className="h-5 w-5" />,
  },
  {
    id: 2,
    name: 'Golden Combination',
    description: 'Selecting best linker-payload combination',
    agent: 'Alchemist',
    icon: <FlaskConical className="h-5 w-5" />,
  },
  {
    id: 3,
    name: 'Property Calculation',
    description: 'Computing molecular properties',
    agent: 'Coder',
    icon: <Calculator className="h-5 w-5" />,
  },
  {
    id: 4,
    name: 'Physical Validation',
    description: 'Verifying structural feasibility',
    agent: 'Auditor',
    icon: <Shield className="h-5 w-5" />,
  },
  {
    id: 5,
    name: 'Virtual Trial',
    description: 'Running clinical simulation',
    agent: 'Clinical',
    icon: <Activity className="h-5 w-5" />,
  },
];

const statusColors = {
  pending: 'text-slate-500 bg-slate-800',
  running: 'text-violet-400 bg-violet-900/50',
  completed: 'text-emerald-400 bg-emerald-900/50',
  error: 'text-red-400 bg-red-900/50',
};

const statusIcons = {
  pending: <Circle className="h-5 w-5" />,
  running: <Loader2 className="h-5 w-5 animate-spin" />,
  completed: <CheckCircle2 className="h-5 w-5" />,
  error: <AlertCircle className="h-5 w-5" />,
};

export function NavigatorPipeline({
  currentStep,
  steps,
  diseaseName,
  agentLogs = [],
  onComplete,
}: NavigatorPipelineProps) {
  const [progress, setProgress] = useState(0);

  // Update progress based on current step
  useEffect(() => {
    const completedSteps = steps.filter((s) => s.status === 'completed').length;
    setProgress((completedSteps / DEFAULT_STEPS.length) * 100);

    // Notify completion
    if (completedSteps === DEFAULT_STEPS.length) {
      onComplete?.();
    }
  }, [steps, onComplete]);

  // Combine default steps with current status
  const pipelineSteps: PipelineStep[] = DEFAULT_STEPS.map((step, index) => ({
    ...step,
    status: steps[index]?.status || 'pending',
    message: steps[index]?.message,
  }));

  return (
    <Card className="bg-slate-900/50 border-slate-700 overflow-hidden">
      <CardHeader className="border-b border-slate-700/50">
        <div className="flex items-center justify-between">
          <div>
            <CardTitle className="text-white flex items-center gap-2">
              Live Pipeline
              {currentStep > 0 && currentStep <= 5 && (
                <Badge variant="outline" className="text-violet-400 border-violet-400/50 animate-pulse">
                  Running
                </Badge>
              )}
            </CardTitle>
            <p className="text-sm text-slate-500 mt-1">
              Designing ADC for: <span className="text-violet-400">{diseaseName}</span>
            </p>
          </div>
          <div className="text-right">
            <p className="text-sm text-slate-500">Progress</p>
            <p className="text-2xl font-bold text-white">{Math.round(progress)}%</p>
          </div>
        </div>
        <Progress value={progress} className="mt-4 h-2" />
      </CardHeader>

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

      {/* FIXED: 실시간 에이전트 로그 패널 */}
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
