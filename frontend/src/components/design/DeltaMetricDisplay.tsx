/**
 * Delta Metric Display - 전후 비교 화살표 지표
 * Lead Optimization 최적화 효과 시각화
 */
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import { Progress } from '@/components/ui/progress'
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger
} from '@/components/ui/tooltip'
import { motion } from 'framer-motion'
import {
  TrendingUp,
  TrendingDown,
  Minus,
  ArrowRight,
  Activity,
  AlertTriangle,
  CheckCircle2,
  Info
} from 'lucide-react'

export interface MetricDelta {
  name: string
  label: string
  before: number
  after: number
  unit?: string
  idealDirection: 'increase' | 'decrease' | 'neutral'
  threshold?: {
    good: number
    warning: number
  }
  description?: string
}

interface DeltaMetricDisplayProps {
  metrics: MetricDelta[]
  title?: string
  showSummary?: boolean
  className?: string
}

export function DeltaMetricDisplay({
  metrics,
  title = 'Optimization Results',
  showSummary = true,
  className
}: DeltaMetricDisplayProps) {
  // Calculate improvement score
  const improvements = metrics.filter(m => {
    const delta = m.after - m.before
    return (m.idealDirection === 'increase' && delta > 0) ||
           (m.idealDirection === 'decrease' && delta < 0)
  })
  const improvementRate = (improvements.length / metrics.length) * 100

  return (
    <Card className={className}>
      <CardHeader className="pb-2">
        <div className="flex items-center justify-between">
          <CardTitle className="text-sm font-medium flex items-center gap-2">
            <Activity className="w-4 h-4 text-blue-500" />
            {title}
          </CardTitle>
          {showSummary && (
            <Badge
              variant="outline"
              className={`text-xs ${
                improvementRate >= 70
                  ? 'text-green-600 border-green-200 bg-green-50'
                  : improvementRate >= 50
                  ? 'text-amber-600 border-amber-200 bg-amber-50'
                  : 'text-red-600 border-red-200 bg-red-50'
              }`}
            >
              {improvements.length}/{metrics.length} Improved
            </Badge>
          )}
        </div>
      </CardHeader>
      <CardContent className="space-y-3">
        {metrics.map((metric, index) => (
          <MetricDeltaRow key={metric.name} metric={metric} index={index} />
        ))}
      </CardContent>
    </Card>
  )
}

function MetricDeltaRow({
  metric,
  index
}: {
  metric: MetricDelta
  index: number
}) {
  const delta = metric.after - metric.before
  const deltaPercent = metric.before !== 0
    ? ((delta / Math.abs(metric.before)) * 100)
    : 0

  const isImproved =
    (metric.idealDirection === 'increase' && delta > 0) ||
    (metric.idealDirection === 'decrease' && delta < 0)
  const isWorsened =
    (metric.idealDirection === 'increase' && delta < 0) ||
    (metric.idealDirection === 'decrease' && delta > 0)

  const getDeltaColor = () => {
    if (isImproved) return 'text-green-600'
    if (isWorsened) return 'text-red-600'
    return 'text-gray-500'
  }

  const getDeltaIcon = () => {
    if (delta > 0) return <TrendingUp className="w-4 h-4" />
    if (delta < 0) return <TrendingDown className="w-4 h-4" />
    return <Minus className="w-4 h-4" />
  }

  const getStatusIcon = () => {
    if (isImproved) return <CheckCircle2 className="w-4 h-4 text-green-500" />
    if (isWorsened) return <AlertTriangle className="w-4 h-4 text-red-500" />
    return <Minus className="w-4 h-4 text-gray-400" />
  }

  // Progress bar calculation (normalized to 0-100)
  const normalizeValue = (val: number) => {
    if (metric.threshold) {
      const max = metric.threshold.warning * 1.5
      return Math.min(100, (val / max) * 100)
    }
    return Math.min(100, val)
  }

  return (
    <motion.div
      initial={{ opacity: 0, x: -10 }}
      animate={{ opacity: 1, x: 0 }}
      transition={{ delay: index * 0.05 }}
      className="p-3 rounded-lg border bg-slate-50/50"
    >
      <div className="flex items-center justify-between mb-2">
        <div className="flex items-center gap-2">
          <TooltipProvider>
            <Tooltip>
              <TooltipTrigger asChild>
                <span className="text-xs font-medium flex items-center gap-1 cursor-help">
                  {metric.label}
                  {metric.description && (
                    <Info className="w-3 h-3 text-gray-400" />
                  )}
                </span>
              </TooltipTrigger>
              {metric.description && (
                <TooltipContent>
                  <p className="text-xs max-w-[200px]">{metric.description}</p>
                </TooltipContent>
              )}
            </Tooltip>
          </TooltipProvider>
        </div>
        {getStatusIcon()}
      </div>

      {/* Before -> After with Delta */}
      <div className="flex items-center gap-3">
        {/* Before Value */}
        <div className="text-center min-w-[60px]">
          <p className="text-[10px] text-gray-500 uppercase">Before</p>
          <p className="text-sm font-mono font-medium text-gray-600">
            {metric.before.toFixed(2)}
            {metric.unit && <span className="text-[10px] ml-0.5">{metric.unit}</span>}
          </p>
        </div>

        {/* Arrow with Delta */}
        <div className="flex-1 flex items-center justify-center">
          <div className={`flex items-center gap-1 px-2 py-1 rounded-full ${
            isImproved ? 'bg-green-100' :
            isWorsened ? 'bg-red-100' : 'bg-gray-100'
          }`}>
            <span className={`text-xs font-medium ${getDeltaColor()}`}>
              {delta >= 0 ? '+' : ''}{delta.toFixed(2)}
            </span>
            {getDeltaIcon()}
          </div>
        </div>

        {/* After Value */}
        <div className="text-center min-w-[60px]">
          <p className="text-[10px] text-gray-500 uppercase">After</p>
          <p className={`text-sm font-mono font-semibold ${
            isImproved ? 'text-green-600' :
            isWorsened ? 'text-red-600' : 'text-gray-800'
          }`}>
            {metric.after.toFixed(2)}
            {metric.unit && <span className="text-[10px] ml-0.5">{metric.unit}</span>}
          </p>
        </div>
      </div>

      {/* Progress Bars (optional visual) */}
      <div className="mt-2 space-y-1">
        <div className="flex items-center gap-2">
          <span className="text-[9px] text-gray-400 w-8">Before</span>
          <Progress
            value={normalizeValue(metric.before)}
            className="h-1.5 flex-1"
          />
        </div>
        <div className="flex items-center gap-2">
          <span className="text-[9px] text-gray-400 w-8">After</span>
          <Progress
            value={normalizeValue(metric.after)}
            className={`h-1.5 flex-1 ${
              isImproved ? '[&>div]:bg-green-500' :
              isWorsened ? '[&>div]:bg-red-500' : ''
            }`}
          />
        </div>
      </div>

      {/* Percentage Change */}
      <div className="mt-2 flex justify-end">
        <Badge
          variant="secondary"
          className={`text-[10px] ${getDeltaColor()}`}
        >
          {deltaPercent >= 0 ? '+' : ''}{deltaPercent.toFixed(1)}%
        </Badge>
      </div>
    </motion.div>
  )
}

/**
 * Compact inline delta display
 */
interface CompactDeltaProps {
  before: number
  after: number
  unit?: string
  idealDirection?: 'increase' | 'decrease' | 'neutral'
}

export function CompactDelta({
  before,
  after,
  unit = '',
  idealDirection = 'neutral'
}: CompactDeltaProps) {
  const delta = after - before
  const isImproved =
    (idealDirection === 'increase' && delta > 0) ||
    (idealDirection === 'decrease' && delta < 0)
  const isWorsened =
    (idealDirection === 'increase' && delta < 0) ||
    (idealDirection === 'decrease' && delta > 0)

  return (
    <span className="inline-flex items-center gap-1 text-xs">
      <span className="text-gray-500">{before.toFixed(1)}</span>
      <ArrowRight className="w-3 h-3 text-gray-400" />
      <span className={
        isImproved ? 'text-green-600 font-medium' :
        isWorsened ? 'text-red-600 font-medium' : 'text-gray-700'
      }>
        {after.toFixed(1)}{unit}
      </span>
      <span className={`text-[10px] ${
        delta >= 0 ? 'text-green-600' : 'text-red-600'
      }`}>
        ({delta >= 0 ? '+' : ''}{delta.toFixed(1)})
      </span>
    </span>
  )
}

/**
 * Create common ADC optimization metrics
 */
export function createOptimizationMetrics(data: {
  toxicityBefore?: number
  toxicityAfter?: number
  stabilityBefore?: number
  stabilityAfter?: number
  saScoreBefore?: number
  saScoreAfter?: number
  logPBefore?: number
  logPAfter?: number
  mwBefore?: number
  mwAfter?: number
}): MetricDelta[] {
  const metrics: MetricDelta[] = []

  if (data.toxicityBefore !== undefined && data.toxicityAfter !== undefined) {
    metrics.push({
      name: 'toxicity',
      label: 'Toxicity Score',
      before: data.toxicityBefore,
      after: data.toxicityAfter,
      idealDirection: 'decrease',
      description: 'Lower is better. Predicted off-target toxicity.'
    })
  }

  if (data.stabilityBefore !== undefined && data.stabilityAfter !== undefined) {
    metrics.push({
      name: 'stability',
      label: 'Metabolic Stability',
      before: data.stabilityBefore,
      after: data.stabilityAfter,
      unit: '%',
      idealDirection: 'increase',
      description: 'Higher is better. % remaining after 1hr in HLM.'
    })
  }

  if (data.saScoreBefore !== undefined && data.saScoreAfter !== undefined) {
    metrics.push({
      name: 'saScore',
      label: 'SA Score',
      before: data.saScoreBefore,
      after: data.saScoreAfter,
      idealDirection: 'decrease',
      threshold: { good: 3, warning: 5 },
      description: 'Synthesis Accessibility. 1=easy, 10=hard.'
    })
  }

  if (data.logPBefore !== undefined && data.logPAfter !== undefined) {
    metrics.push({
      name: 'logP',
      label: 'LogP',
      before: data.logPBefore,
      after: data.logPAfter,
      idealDirection: 'neutral',
      description: 'Lipophilicity. Ideal range: 2-5.'
    })
  }

  if (data.mwBefore !== undefined && data.mwAfter !== undefined) {
    metrics.push({
      name: 'mw',
      label: 'Molecular Weight',
      before: data.mwBefore,
      after: data.mwAfter,
      unit: 'Da',
      idealDirection: 'neutral',
      description: 'Ideal for payload: 300-500 Da.'
    })
  }

  return metrics
}
