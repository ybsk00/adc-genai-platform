/**
 * Auditor Feedback - 반려 및 피드백 UI
 * Redesign Loop 상태 시각화
 */
import { useState } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import { Button } from '@/components/ui/button'
import { Progress } from '@/components/ui/progress'
import {
  Popover,
  PopoverContent,
  PopoverTrigger
} from '@/components/ui/popover'
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger
} from '@/components/ui/tooltip'
import { motion, AnimatePresence } from 'framer-motion'
import {
  Shield,
  ShieldAlert,
  ShieldCheck,
  AlertTriangle,
  CheckCircle2,
  XCircle,
  RefreshCcw,
  Info,
  ChevronDown,
  ArrowRight,
  Target,
  Beaker
} from 'lucide-react'

export interface AuditCheckItem {
  id: string
  name: string
  category: 'lipinski' | 'pains' | 'stability' | 'toxicity' | 'dar' | 'other'
  status: 'pass' | 'fail' | 'warning' | 'pending'
  currentValue?: number | string
  targetValue?: number | string
  unit?: string
  reasoning?: string
}

export interface RedesignRequest {
  requestedBy: 'auditor'
  reason: string
  failedChecks: string[]
  suggestions: string[]
  targetAgent: 'alchemist' | 'coder'
  iteration: number
  maxIterations: number
  timestamp: string
}

interface AuditorFeedbackProps {
  checkItems: AuditCheckItem[]
  redesignRequest?: RedesignRequest | null
  isInRedesignLoop?: boolean
  onRedesignConfirm?: () => void
  className?: string
}

export function AuditorFeedback({
  checkItems,
  redesignRequest,
  isInRedesignLoop = false,
  onRedesignConfirm,
  className
}: AuditorFeedbackProps) {
  const [expandedCategory, setExpandedCategory] = useState<string | null>(null)

  // Group items by category
  const groupedItems = checkItems.reduce((acc, item) => {
    if (!acc[item.category]) acc[item.category] = []
    acc[item.category].push(item)
    return acc
  }, {} as Record<string, AuditCheckItem[]>)

  // Calculate overall status
  const passCount = checkItems.filter(i => i.status === 'pass').length
  const failCount = checkItems.filter(i => i.status === 'fail').length
  const warningCount = checkItems.filter(i => i.status === 'warning').length
  const overallPass = failCount === 0

  const categoryLabels: Record<string, string> = {
    lipinski: 'Lipinski Rule of 5',
    pains: 'PAINS Filter',
    stability: 'Metabolic Stability',
    toxicity: 'Toxicity Assessment',
    dar: 'DAR Compliance',
    other: 'Other Checks'
  }

  return (
    <Card className={className}>
      <CardHeader className="pb-2">
        <div className="flex items-center justify-between">
          <CardTitle className="text-sm font-medium flex items-center gap-2">
            <Shield className="w-4 h-4 text-red-500" />
            Audit Checklist
          </CardTitle>
          <Badge
            variant="outline"
            className={`text-xs ${
              overallPass
                ? 'text-green-600 border-green-200 bg-green-50'
                : 'text-red-600 border-red-200 bg-red-50'
            }`}
          >
            {overallPass ? (
              <>
                <ShieldCheck className="w-3 h-3 mr-1" />
                All Passed
              </>
            ) : (
              <>
                <ShieldAlert className="w-3 h-3 mr-1" />
                {failCount} Failed
              </>
            )}
          </Badge>
        </div>

        {/* Summary Bar */}
        <div className="flex items-center gap-2 mt-2">
          <Progress
            value={(passCount / checkItems.length) * 100}
            className="flex-1 h-2"
          />
          <span className="text-[10px] text-gray-500 whitespace-nowrap">
            {passCount}/{checkItems.length}
          </span>
        </div>

        {/* Status Pills */}
        <div className="flex gap-2 mt-2">
          <Badge variant="secondary" className="text-[10px] bg-green-50 text-green-600">
            <CheckCircle2 className="w-3 h-3 mr-1" />
            {passCount} Pass
          </Badge>
          <Badge variant="secondary" className="text-[10px] bg-red-50 text-red-600">
            <XCircle className="w-3 h-3 mr-1" />
            {failCount} Fail
          </Badge>
          <Badge variant="secondary" className="text-[10px] bg-amber-50 text-amber-600">
            <AlertTriangle className="w-3 h-3 mr-1" />
            {warningCount} Warning
          </Badge>
        </div>
      </CardHeader>

      <CardContent className="space-y-3">
        {/* Redesign Loop Alert */}
        <AnimatePresence>
          {isInRedesignLoop && redesignRequest && (
            <RedesignLoopAlert
              request={redesignRequest}
              onConfirm={onRedesignConfirm}
            />
          )}
        </AnimatePresence>

        {/* Category Groups */}
        {Object.entries(groupedItems).map(([category, items]) => (
          <CategoryGroup
            key={category}
            category={category}
            label={categoryLabels[category] || category}
            items={items}
            isExpanded={expandedCategory === category}
            onToggle={() => setExpandedCategory(
              expandedCategory === category ? null : category
            )}
          />
        ))}
      </CardContent>
    </Card>
  )
}

// Redesign Loop Alert Component
function RedesignLoopAlert({
  request,
  onConfirm
}: {
  request: RedesignRequest
  onConfirm?: () => void
}) {
  return (
    <motion.div
      initial={{ opacity: 0, scale: 0.95 }}
      animate={{ opacity: 1, scale: 1 }}
      exit={{ opacity: 0, scale: 0.95 }}
      className="p-3 rounded-lg bg-amber-50 border-2 border-amber-300"
    >
      <div className="flex items-start gap-2">
        <motion.div
          animate={{ rotate: 360 }}
          transition={{ repeat: Infinity, duration: 2, ease: 'linear' }}
        >
          <RefreshCcw className="w-5 h-5 text-amber-600" />
        </motion.div>
        <div className="flex-1">
          <h4 className="text-sm font-semibold text-amber-800 flex items-center gap-2">
            Redesign Loop Active
            <Badge variant="outline" className="text-[10px] text-amber-600">
              Iteration {request.iteration}/{request.maxIterations}
            </Badge>
          </h4>

          <p className="text-xs text-amber-700 mt-1">
            {request.reason}
          </p>

          {/* Failed Checks */}
          <div className="mt-2 flex flex-wrap gap-1">
            {request.failedChecks.map((check, idx) => (
              <Badge
                key={idx}
                variant="destructive"
                className="text-[10px]"
              >
                {check}
              </Badge>
            ))}
          </div>

          {/* Target Agent */}
          <div className="mt-2 flex items-center gap-2 text-xs text-amber-700">
            <span>Requesting fix from:</span>
            <Badge variant="outline" className="text-amber-600">
              {request.targetAgent === 'alchemist' ? (
                <>
                  <Beaker className="w-3 h-3 mr-1" />
                  Alchemist
                </>
              ) : (
                <>
                  <Target className="w-3 h-3 mr-1" />
                  Coder
                </>
              )}
            </Badge>
          </div>

          {/* Suggestions */}
          {request.suggestions.length > 0 && (
            <div className="mt-2 p-2 bg-white/50 rounded border border-amber-200">
              <p className="text-[10px] text-amber-600 font-medium mb-1">
                Suggestions:
              </p>
              <ul className="text-[11px] text-amber-800 space-y-0.5">
                {request.suggestions.map((s, idx) => (
                  <li key={idx} className="flex items-start gap-1">
                    <ArrowRight className="w-3 h-3 mt-0.5 shrink-0" />
                    {s}
                  </li>
                ))}
              </ul>
            </div>
          )}
        </div>
      </div>
    </motion.div>
  )
}

// Category Group Component
function CategoryGroup({
  category,
  label,
  items,
  isExpanded,
  onToggle
}: {
  category: string
  label: string
  items: AuditCheckItem[]
  isExpanded: boolean
  onToggle: () => void
}) {
  const passCount = items.filter(i => i.status === 'pass').length
  const allPass = passCount === items.length

  return (
    <div className="border rounded-lg overflow-hidden">
      <button
        className={`w-full flex items-center justify-between p-2 text-left hover:bg-slate-50 ${
          allPass ? 'bg-green-50/30' : 'bg-red-50/30'
        }`}
        onClick={onToggle}
      >
        <div className="flex items-center gap-2">
          {allPass ? (
            <CheckCircle2 className="w-4 h-4 text-green-500" />
          ) : (
            <XCircle className="w-4 h-4 text-red-500" />
          )}
          <span className="text-xs font-medium">{label}</span>
        </div>
        <div className="flex items-center gap-2">
          <Badge variant="secondary" className="text-[10px]">
            {passCount}/{items.length}
          </Badge>
          <ChevronDown
            className={`w-4 h-4 text-gray-400 transition-transform ${
              isExpanded ? 'rotate-180' : ''
            }`}
          />
        </div>
      </button>

      <AnimatePresence>
        {isExpanded && (
          <motion.div
            initial={{ height: 0, opacity: 0 }}
            animate={{ height: 'auto', opacity: 1 }}
            exit={{ height: 0, opacity: 0 }}
            className="border-t"
          >
            <div className="p-2 space-y-1">
              {items.map((item) => (
                <CheckItemRow key={item.id} item={item} />
              ))}
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  )
}

// Check Item Row Component
function CheckItemRow({ item }: { item: AuditCheckItem }) {
  const getStatusIcon = () => {
    switch (item.status) {
      case 'pass':
        return <CheckCircle2 className="w-3 h-3 text-green-500" />
      case 'fail':
        return <XCircle className="w-3 h-3 text-red-500" />
      case 'warning':
        return <AlertTriangle className="w-3 h-3 text-amber-500" />
      default:
        return <RefreshCcw className="w-3 h-3 text-gray-400 animate-spin" />
    }
  }

  return (
    <div className={`flex items-center justify-between p-2 rounded ${
      item.status === 'pass' ? 'bg-green-50' :
      item.status === 'fail' ? 'bg-red-50' :
      item.status === 'warning' ? 'bg-amber-50' : 'bg-gray-50'
    }`}>
      <div className="flex items-center gap-2">
        {getStatusIcon()}
        <span className="text-xs">{item.name}</span>

        {/* Reasoning Popover */}
        {item.reasoning && (
          <Popover>
            <PopoverTrigger asChild>
              <Button variant="ghost" size="icon" className="h-5 w-5">
                <Info className="w-3 h-3 text-gray-400" />
              </Button>
            </PopoverTrigger>
            <PopoverContent className="w-64 text-xs">
              <p className="font-medium mb-1">Auditor Reasoning:</p>
              <p className="text-gray-600">{item.reasoning}</p>
            </PopoverContent>
          </Popover>
        )}
      </div>

      {/* Value Display */}
      {item.currentValue !== undefined && (
        <div className="text-right">
          <span className={`text-xs font-mono ${
            item.status === 'pass' ? 'text-green-600' :
            item.status === 'fail' ? 'text-red-600' : 'text-amber-600'
          }`}>
            {item.currentValue}
            {item.unit && <span className="text-[10px] ml-0.5">{item.unit}</span>}
          </span>
          {item.targetValue !== undefined && (
            <span className="text-[10px] text-gray-400 ml-1">
              (target: {item.targetValue})
            </span>
          )}
        </div>
      )}
    </div>
  )
}

/**
 * Compact audit status badge
 */
export function CompactAuditBadge({
  passCount,
  totalCount,
  onClick
}: {
  passCount: number
  totalCount: number
  onClick?: () => void
}) {
  const allPass = passCount === totalCount

  return (
    <Badge
      variant="outline"
      className={`cursor-pointer ${
        allPass
          ? 'text-green-600 border-green-200 hover:bg-green-50'
          : 'text-red-600 border-red-200 hover:bg-red-50'
      }`}
      onClick={onClick}
    >
      {allPass ? (
        <ShieldCheck className="w-3 h-3 mr-1" />
      ) : (
        <ShieldAlert className="w-3 h-3 mr-1" />
      )}
      {passCount}/{totalCount}
    </Badge>
  )
}
