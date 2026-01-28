/**
 * Success Radar Chart - 설계 지표 거미줄 차트
 * MW, LogP, DAR 일치도, Golden Set 유사도 시각화
 */
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import {
  RadarChart,
  PolarGrid,
  PolarAngleAxis,
  PolarRadiusAxis,
  Radar,
  ResponsiveContainer,
  Legend,
  Tooltip
} from 'recharts'
import { Target, TrendingUp } from 'lucide-react'

export interface RadarMetric {
  subject: string
  value: number       // 실제 값 (0-100 스케일)
  target?: number     // 목표 값 (0-100 스케일)
  fullMark: number    // 최대값 (보통 100)
}

interface SuccessRadarChartProps {
  metrics: RadarMetric[]
  title?: string
  showTarget?: boolean
  overallScore?: number
  className?: string
}

// 기본 메트릭 템플릿
export function createDefaultMetrics(data: {
  mw?: number          // Molecular Weight (ideal: 300-500)
  logP?: number        // LogP (ideal: 2-5)
  darMatch?: number    // DAR 일치도 (0-100%)
  goldenSetSim?: number // Golden Set 유사도 (0-100%)
  tpsa?: number        // Topological Polar Surface Area
  hbd?: number         // Hydrogen Bond Donors (ideal: < 5)
}): RadarMetric[] {
  const normalize = (val: number, min: number, max: number): number => {
    if (val < min) return ((val / min) * 50)
    if (val > max) return Math.max(0, 100 - ((val - max) / max) * 50)
    return 50 + ((val - min) / (max - min)) * 50
  }

  return [
    {
      subject: 'MW',
      value: data.mw ? normalize(data.mw, 300, 500) : 0,
      target: 75,
      fullMark: 100
    },
    {
      subject: 'LogP',
      value: data.logP ? normalize(data.logP, 2, 5) : 0,
      target: 80,
      fullMark: 100
    },
    {
      subject: 'DAR Match',
      value: data.darMatch ?? 0,
      target: 90,
      fullMark: 100
    },
    {
      subject: 'Golden Set',
      value: data.goldenSetSim ?? 0,
      target: 70,
      fullMark: 100
    },
    {
      subject: 'TPSA',
      value: data.tpsa ? normalize(data.tpsa, 40, 140) : 0,
      target: 75,
      fullMark: 100
    },
    {
      subject: 'H-Bond',
      value: data.hbd ? Math.min(100, (1 - data.hbd / 10) * 100) : 0,
      target: 80,
      fullMark: 100
    }
  ]
}

export function SuccessRadarChart({
  metrics,
  title = 'Success Metrics',
  showTarget = true,
  overallScore,
  className
}: SuccessRadarChartProps) {
  // Calculate overall score if not provided
  const calculatedScore = overallScore ?? (
    metrics.reduce((sum, m) => sum + m.value, 0) / metrics.length
  )

  const getScoreColor = (score: number) => {
    if (score >= 80) return 'text-green-600 bg-green-50 border-green-200'
    if (score >= 60) return 'text-amber-600 bg-amber-50 border-amber-200'
    return 'text-red-600 bg-red-50 border-red-200'
  }

  const getScoreLabel = (score: number) => {
    if (score >= 80) return 'Excellent'
    if (score >= 60) return 'Good'
    if (score >= 40) return 'Fair'
    return 'Needs Improvement'
  }

  return (
    <Card className={className}>
      <CardHeader className="pb-2">
        <div className="flex items-center justify-between">
          <CardTitle className="text-sm font-medium flex items-center gap-2">
            <Target className="w-4 h-4 text-blue-500" />
            {title}
          </CardTitle>
          <Badge
            variant="outline"
            className={`text-xs ${getScoreColor(calculatedScore)}`}
          >
            <TrendingUp className="w-3 h-3 mr-1" />
            {calculatedScore.toFixed(0)}% - {getScoreLabel(calculatedScore)}
          </Badge>
        </div>
      </CardHeader>
      <CardContent>
        <div className="h-[280px] w-full">
          <ResponsiveContainer width="100%" height="100%">
            <RadarChart cx="50%" cy="50%" outerRadius="70%" data={metrics}>
              <PolarGrid strokeDasharray="3 3" />
              <PolarAngleAxis
                dataKey="subject"
                tick={{ fontSize: 11, fill: '#6b7280' }}
              />
              <PolarRadiusAxis
                angle={30}
                domain={[0, 100]}
                tick={{ fontSize: 10 }}
                tickCount={5}
              />

              {/* Target area (optional) */}
              {showTarget && (
                <Radar
                  name="Target"
                  dataKey="target"
                  stroke="#94a3b8"
                  fill="#94a3b8"
                  fillOpacity={0.1}
                  strokeDasharray="5 5"
                />
              )}

              {/* Actual values */}
              <Radar
                name="Current"
                dataKey="value"
                stroke="#3b82f6"
                fill="#3b82f6"
                fillOpacity={0.4}
                strokeWidth={2}
              />

              <Tooltip
                content={({ payload }) => {
                  if (!payload?.length) return null
                  const data = payload[0]?.payload
                  return (
                    <div className="bg-white border rounded-lg shadow-lg p-2 text-xs">
                      <p className="font-semibold">{data?.subject}</p>
                      <p className="text-blue-600">Value: {data?.value?.toFixed(1)}%</p>
                      {data?.target && (
                        <p className="text-gray-500">Target: {data?.target}%</p>
                      )}
                    </div>
                  )
                }}
              />

              <Legend
                wrapperStyle={{ fontSize: '12px' }}
                iconType="circle"
              />
            </RadarChart>
          </ResponsiveContainer>
        </div>

        {/* Metric breakdown */}
        <div className="grid grid-cols-3 gap-2 mt-2 pt-2 border-t">
          {metrics.slice(0, 6).map((metric) => (
            <div
              key={metric.subject}
              className="text-center p-1 rounded bg-slate-50"
            >
              <p className="text-[10px] text-gray-500 uppercase truncate">
                {metric.subject}
              </p>
              <p className={`text-xs font-semibold ${
                metric.value >= 70 ? 'text-green-600' :
                metric.value >= 50 ? 'text-amber-600' : 'text-red-600'
              }`}>
                {metric.value.toFixed(0)}%
              </p>
            </div>
          ))}
        </div>
      </CardContent>
    </Card>
  )
}
