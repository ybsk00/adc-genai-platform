/**
 * Design Progress Timeline
 * 설계 워크플로우 진행 상황 타임라인
 */
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Progress } from '@/components/ui/progress'
import { motion } from 'framer-motion'
import {
  CheckCircle2,
  Circle,
  Loader2,
  BookOpen,
  Beaker,
  Code2,
  HeartPulse,
  Shield
} from 'lucide-react'

interface TimelineStep {
  id: string
  name: string
  description: string
  icon: React.ReactNode
}

const WORKFLOW_STEPS: TimelineStep[] = [
  {
    id: 'librarian',
    name: 'Knowledge Search',
    description: 'Searching literature and Golden Set references',
    icon: <BookOpen className="w-5 h-5" />
  },
  {
    id: 'alchemist',
    name: 'Structure Design',
    description: 'Generating candidate SMILES structures',
    icon: <Beaker className="w-5 h-5" />
  },
  {
    id: 'coder',
    name: 'Property Calculation',
    description: 'Computing molecular properties',
    icon: <Code2 className="w-5 h-5" />
  },
  {
    id: 'auditor',
    name: 'Validation',
    description: 'Chemistry validation and constraint check',
    icon: <Shield className="w-5 h-5" />
  }
]

interface DesignProgressTimelineProps {
  currentStep: number
  currentAgent: string | null
  status: string | null
  healingAttempts?: number
}

export function DesignProgressTimeline({
  currentStep,
  currentAgent,
  status,
  healingAttempts = 0
}: DesignProgressTimelineProps) {
  const progressPercent = Math.min(100, (currentStep / WORKFLOW_STEPS.length) * 100)

  const getStepStatus = (stepIndex: number) => {
    if (status === 'completed' || status === 'manual_review') {
      return 'completed'
    }
    if (stepIndex < currentStep) {
      return 'completed'
    }
    if (stepIndex === currentStep) {
      return 'active'
    }
    return 'pending'
  }

  const isHealing = currentAgent === 'healer'

  return (
    <Card>
      <CardHeader className="pb-2">
        <CardTitle className="text-sm font-medium flex items-center justify-between">
          <span>Design Progress</span>
          <span className="text-xs font-normal text-gray-500">
            Step {currentStep} / {WORKFLOW_STEPS.length}
          </span>
        </CardTitle>
      </CardHeader>
      <CardContent className="space-y-4">
        {/* Progress bar */}
        <div className="space-y-1">
          <Progress value={progressPercent} className="h-2" />
          <div className="flex justify-between text-xs text-gray-500">
            <span>
              {status === 'completed'
                ? 'Completed'
                : status === 'failed'
                ? 'Failed'
                : status === 'manual_review'
                ? 'Manual Review Required'
                : 'In Progress'}
            </span>
            <span>{Math.round(progressPercent)}%</span>
          </div>
        </div>

        {/* Healing indicator */}
        {isHealing && (
          <motion.div
            initial={{ opacity: 0, y: -10 }}
            animate={{ opacity: 1, y: 0 }}
            className="flex items-center gap-2 p-2 bg-amber-50 border border-amber-200 rounded-lg text-sm"
          >
            <HeartPulse className="w-4 h-4 text-amber-500" />
            <span className="text-amber-700">
              Self-healing in progress (Attempt {healingAttempts}/3)
            </span>
          </motion.div>
        )}

        {/* Timeline steps */}
        <div className="relative">
          {/* Vertical line */}
          <div className="absolute left-[15px] top-0 bottom-0 w-0.5 bg-gray-200" />

          <div className="space-y-4">
            {WORKFLOW_STEPS.map((step, index) => {
              const stepStatus = getStepStatus(index + 1)
              const isActive = stepStatus === 'active' && !isHealing

              return (
                <motion.div
                  key={step.id}
                  initial={{ opacity: 0, x: -10 }}
                  animate={{ opacity: 1, x: 0 }}
                  transition={{ delay: index * 0.1 }}
                  className={`relative flex items-start gap-3 ${
                    stepStatus === 'pending' ? 'opacity-50' : ''
                  }`}
                >
                  {/* Step indicator */}
                  <div
                    className={`relative z-10 flex items-center justify-center w-8 h-8 rounded-full border-2 ${
                      stepStatus === 'completed'
                        ? 'bg-green-500 border-green-500 text-white'
                        : isActive
                        ? 'bg-blue-500 border-blue-500 text-white'
                        : 'bg-white border-gray-300 text-gray-400'
                    }`}
                  >
                    {stepStatus === 'completed' ? (
                      <CheckCircle2 className="w-4 h-4" />
                    ) : isActive ? (
                      <Loader2 className="w-4 h-4 animate-spin" />
                    ) : (
                      step.icon
                    )}
                  </div>

                  {/* Step content */}
                  <div className="flex-1 min-w-0 pt-0.5">
                    <p
                      className={`text-sm font-medium ${
                        stepStatus === 'completed'
                          ? 'text-green-700'
                          : isActive
                          ? 'text-blue-700'
                          : 'text-gray-500'
                      }`}
                    >
                      {step.name}
                    </p>
                    <p className="text-xs text-gray-500 mt-0.5">
                      {step.description}
                    </p>
                  </div>
                </motion.div>
              )
            })}
          </div>
        </div>
      </CardContent>
    </Card>
  )
}
