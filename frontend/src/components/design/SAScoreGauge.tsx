/**
 * SA Score Gauge - 합성 용이성 게이지
 * Synthesis Accessibility Score 시각화 (1~10점)
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
  FlaskConical,
  Thermometer,
  AlertTriangle,
  CheckCircle2,
  Info,
  TrendingDown,
  TrendingUp
} from 'lucide-react'

interface SAScoreGaugeProps {
  score: number  // 1-10 scale
  previousScore?: number
  details?: SAScoreDetails
  className?: string
}

interface SAScoreDetails {
  fragmentScore?: number
  ringScore?: number
  stereoScore?: number
  sizeScore?: number
  bridgeScore?: number
}

export function SAScoreGauge({
  score,
  previousScore,
  details,
  className
}: SAScoreGaugeProps) {
  // Normalize score to 0-100 for progress bar
  const normalizedScore = Math.max(0, Math.min(100, (10 - score) * 10))

  // Get color based on score
  const getScoreConfig = (s: number) => {
    if (s <= 3) return {
      label: 'Easy',
      color: 'text-green-600',
      bgColor: 'bg-green-50 border-green-200',
      progressColor: 'bg-green-500'
    }
    if (s <= 5) return {
      label: 'Moderate',
      color: 'text-amber-600',
      bgColor: 'bg-amber-50 border-amber-200',
      progressColor: 'bg-amber-500'
    }
    if (s <= 7) return {
      label: 'Difficult',
      color: 'text-orange-600',
      bgColor: 'bg-orange-50 border-orange-200',
      progressColor: 'bg-orange-500'
    }
    return {
      label: 'Very Hard',
      color: 'text-red-600',
      bgColor: 'bg-red-50 border-red-200',
      progressColor: 'bg-red-500'
    }
  }

  const config = getScoreConfig(score)

  // Calculate delta if previous score exists
  const delta = previousScore !== undefined ? previousScore - score : undefined
  const isImproved = delta !== undefined && delta > 0

  return (
    <Card className={className}>
      <CardHeader className="pb-2">
        <div className="flex items-center justify-between">
          <CardTitle className="text-sm font-medium flex items-center gap-2">
            <FlaskConical className="w-4 h-4 text-purple-500" />
            Synthesis Accessibility
          </CardTitle>
          <TooltipProvider>
            <Tooltip>
              <TooltipTrigger asChild>
                <Badge variant="outline" className="text-xs cursor-help">
                  <Info className="w-3 h-3" />
                </Badge>
              </TooltipTrigger>
              <TooltipContent className="max-w-[200px]">
                <p className="text-xs">
                  SA Score (1-10): Lower is easier to synthesize.
                  1 = trivial, 10 = extremely difficult
                </p>
              </TooltipContent>
            </Tooltip>
          </TooltipProvider>
        </div>
      </CardHeader>

      <CardContent className="space-y-4">
        {/* Main Score Display */}
        <div className={`p-4 rounded-lg border ${config.bgColor}`}>
          <div className="flex items-center justify-between">
            <div>
              <motion.div
                key={score}
                initial={{ scale: 1.2, opacity: 0 }}
                animate={{ scale: 1, opacity: 1 }}
                className={`text-4xl font-bold ${config.color}`}
              >
                {score.toFixed(1)}
              </motion.div>
              <p className={`text-xs ${config.color} mt-1`}>
                out of 10
              </p>
            </div>

            <div className="text-right">
              <Badge
                variant="outline"
                className={`text-sm ${config.color} border-current`}
              >
                {config.label}
              </Badge>

              {/* Delta indicator */}
              {delta !== undefined && (
                <div className={`flex items-center justify-end gap-1 mt-2 text-xs ${
                  isImproved ? 'text-green-600' : 'text-red-600'
                }`}>
                  {isImproved ? (
                    <TrendingDown className="w-3 h-3" />
                  ) : (
                    <TrendingUp className="w-3 h-3" />
                  )}
                  <span>
                    {delta > 0 ? '-' : '+'}{Math.abs(delta).toFixed(1)}
                  </span>
                </div>
              )}
            </div>
          </div>

          {/* Visual Gauge */}
          <div className="mt-4">
            <div className="flex justify-between text-[10px] text-gray-500 mb-1">
              <span>Easy (1)</span>
              <span>Hard (10)</span>
            </div>
            <div className="relative h-3 bg-gradient-to-r from-green-200 via-amber-200 to-red-200 rounded-full overflow-hidden">
              <motion.div
                initial={{ left: '50%' }}
                animate={{ left: `${(score / 10) * 100}%` }}
                transition={{ type: 'spring', damping: 15 }}
                className="absolute top-0 -translate-x-1/2 w-4 h-4 -mt-0.5"
              >
                <div className="w-4 h-4 bg-white border-2 border-gray-800 rounded-full shadow" />
              </motion.div>
            </div>
          </div>
        </div>

        {/* Score Breakdown */}
        {details && (
          <div className="space-y-2">
            <p className="text-[10px] text-gray-500 font-medium uppercase">
              Score Breakdown
            </p>
            <div className="space-y-1.5">
              {details.fragmentScore !== undefined && (
                <ScoreDetailRow
                  label="Fragment Complexity"
                  value={details.fragmentScore}
                  maxValue={10}
                />
              )}
              {details.ringScore !== undefined && (
                <ScoreDetailRow
                  label="Ring Systems"
                  value={details.ringScore}
                  maxValue={10}
                />
              )}
              {details.stereoScore !== undefined && (
                <ScoreDetailRow
                  label="Stereochemistry"
                  value={details.stereoScore}
                  maxValue={10}
                />
              )}
              {details.sizeScore !== undefined && (
                <ScoreDetailRow
                  label="Molecular Size"
                  value={details.sizeScore}
                  maxValue={10}
                />
              )}
              {details.bridgeScore !== undefined && (
                <ScoreDetailRow
                  label="Bridged/Spiro"
                  value={details.bridgeScore}
                  maxValue={10}
                />
              )}
            </div>
          </div>
        )}

        {/* Interpretation Guide */}
        <div className="grid grid-cols-4 gap-1 pt-2 border-t">
          <ScoreLegendItem score="1-3" label="Easy" color="bg-green-500" />
          <ScoreLegendItem score="4-5" label="Moderate" color="bg-amber-500" />
          <ScoreLegendItem score="6-7" label="Difficult" color="bg-orange-500" />
          <ScoreLegendItem score="8-10" label="Very Hard" color="bg-red-500" />
        </div>
      </CardContent>
    </Card>
  )
}

function ScoreDetailRow({
  label,
  value,
  maxValue
}: {
  label: string
  value: number
  maxValue: number
}) {
  const percent = (value / maxValue) * 100

  return (
    <div className="flex items-center gap-2">
      <span className="text-[10px] text-gray-600 w-28 truncate">{label}</span>
      <div className="flex-1 h-1.5 bg-gray-100 rounded-full overflow-hidden">
        <motion.div
          initial={{ width: 0 }}
          animate={{ width: `${percent}%` }}
          className={`h-full rounded-full ${
            value <= 3 ? 'bg-green-500' :
            value <= 5 ? 'bg-amber-500' :
            value <= 7 ? 'bg-orange-500' : 'bg-red-500'
          }`}
        />
      </div>
      <span className="text-[10px] font-mono text-gray-500 w-8 text-right">
        {value.toFixed(1)}
      </span>
    </div>
  )
}

function ScoreLegendItem({
  score,
  label,
  color
}: {
  score: string
  label: string
  color: string
}) {
  return (
    <div className="text-center">
      <div className={`w-full h-1.5 ${color} rounded-full mb-1`} />
      <p className="text-[9px] text-gray-500">{score}</p>
      <p className="text-[8px] text-gray-400">{label}</p>
    </div>
  )
}

/**
 * Compact SA Score Badge
 */
export function CompactSAScore({
  score,
  showLabel = true
}: {
  score: number
  showLabel?: boolean
}) {
  const getConfig = (s: number) => {
    if (s <= 3) return { label: 'Easy', color: 'text-green-600 bg-green-50' }
    if (s <= 5) return { label: 'Moderate', color: 'text-amber-600 bg-amber-50' }
    if (s <= 7) return { label: 'Difficult', color: 'text-orange-600 bg-orange-50' }
    return { label: 'Hard', color: 'text-red-600 bg-red-50' }
  }

  const config = getConfig(score)

  return (
    <Badge variant="secondary" className={`text-xs ${config.color}`}>
      <FlaskConical className="w-3 h-3 mr-1" />
      SA: {score.toFixed(1)}
      {showLabel && <span className="ml-1">({config.label})</span>}
    </Badge>
  )
}
