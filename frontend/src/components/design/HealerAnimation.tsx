/**
 * Healer Animation - 자가 치유 라이브 중계
 * 코드 에러 수정 애니메이션 시각화
 */
import { useState, useEffect, useRef } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import { Progress } from '@/components/ui/progress'
import { ScrollArea } from '@/components/ui/scroll-area'
import { motion, AnimatePresence } from 'framer-motion'
import {
  HeartPulse,
  AlertTriangle,
  Wrench,
  CheckCircle2,
  Code2,
  Terminal,
  Zap,
  RefreshCcw
} from 'lucide-react'

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

interface HealerAnimationProps {
  action: HealerAction | null
  isActive: boolean
  className?: string
}

export function HealerAnimation({
  action,
  isActive,
  className
}: HealerAnimationProps) {
  const [displayedCode, setDisplayedCode] = useState<string>('')
  const [isTyping, setIsTyping] = useState(false)
  const codeRef = useRef<HTMLPreElement>(null)

  // Typewriter effect for fixed code
  useEffect(() => {
    if (action?.action === 'healing' && action.fixedCode) {
      setIsTyping(true)
      setDisplayedCode('')

      let index = 0
      const code = action.fixedCode
      const interval = setInterval(() => {
        if (index < code.length) {
          setDisplayedCode(prev => prev + code[index])
          index++
        } else {
          setIsTyping(false)
          clearInterval(interval)
        }
      }, 15) // Fast typing speed

      return () => clearInterval(interval)
    }
  }, [action?.fixedCode, action?.action])

  // Auto-scroll code view
  useEffect(() => {
    if (codeRef.current) {
      codeRef.current.scrollTop = codeRef.current.scrollHeight
    }
  }, [displayedCode])

  if (!isActive && !action) {
    return null
  }

  const getStatusConfig = () => {
    switch (action?.action) {
      case 'detecting':
        return {
          icon: <AlertTriangle className="w-5 h-5" />,
          color: 'text-red-500',
          bgColor: 'bg-red-50 border-red-200',
          label: 'Error Detected',
          progress: 20
        }
      case 'analyzing':
        return {
          icon: <RefreshCcw className="w-5 h-5 animate-spin" />,
          color: 'text-amber-500',
          bgColor: 'bg-amber-50 border-amber-200',
          label: 'Analyzing Error',
          progress: 40
        }
      case 'healing':
        return {
          icon: <Wrench className="w-5 h-5 animate-pulse" />,
          color: 'text-blue-500',
          bgColor: 'bg-blue-50 border-blue-200',
          label: 'Healing Code',
          progress: 70
        }
      case 'healed':
        return {
          icon: <CheckCircle2 className="w-5 h-5" />,
          color: 'text-green-500',
          bgColor: 'bg-green-50 border-green-200',
          label: 'Healed Successfully',
          progress: 100
        }
      case 'failed':
        return {
          icon: <AlertTriangle className="w-5 h-5" />,
          color: 'text-red-500',
          bgColor: 'bg-red-50 border-red-200',
          label: 'Healing Failed',
          progress: 100
        }
      default:
        return {
          icon: <HeartPulse className="w-5 h-5" />,
          color: 'text-gray-500',
          bgColor: 'bg-gray-50 border-gray-200',
          label: 'Monitoring',
          progress: 0
        }
    }
  }

  const config = getStatusConfig()

  return (
    <AnimatePresence>
      <motion.div
        initial={{ opacity: 0, y: -20 }}
        animate={{ opacity: 1, y: 0 }}
        exit={{ opacity: 0, y: -20 }}
        className={className}
      >
        <Card className={`border-2 ${config.bgColor}`}>
          <CardHeader className="pb-2">
            <div className="flex items-center justify-between">
              <CardTitle className="text-sm font-medium flex items-center gap-2">
                <motion.div
                  animate={{ scale: [1, 1.2, 1] }}
                  transition={{ repeat: Infinity, duration: 1.5 }}
                  className={config.color}
                >
                  <HeartPulse className="w-4 h-4" />
                </motion.div>
                Self-Healing System
              </CardTitle>
              <Badge variant="outline" className={`${config.color} text-xs`}>
                {config.icon}
                <span className="ml-1">{config.label}</span>
              </Badge>
            </div>

            {/* Progress Bar */}
            <div className="mt-2 space-y-1">
              <Progress value={config.progress} className="h-2" />
              {action?.attempt && action?.maxAttempts && (
                <p className="text-[10px] text-gray-500 text-right">
                  Attempt {action.attempt}/{action.maxAttempts}
                </p>
              )}
            </div>
          </CardHeader>

          <CardContent className="space-y-3">
            {/* Error Type Badge */}
            {action?.errorType && (
              <div className="flex items-center gap-2">
                <Badge variant="destructive" className="text-xs">
                  <Zap className="w-3 h-3 mr-1" />
                  {action.errorType}
                </Badge>
              </div>
            )}

            {/* Error Message */}
            {action?.errorMessage && action.action !== 'healed' && (
              <motion.div
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                className="p-2 bg-red-100 border border-red-300 rounded-md"
              >
                <p className="text-[10px] text-red-600 font-medium mb-1">
                  <AlertTriangle className="w-3 h-3 inline mr-1" />
                  Error Message
                </p>
                <code className="text-[11px] text-red-800 font-mono block whitespace-pre-wrap break-all">
                  {action.errorMessage.slice(0, 200)}
                  {action.errorMessage.length > 200 && '...'}
                </code>
              </motion.div>
            )}

            {/* Code View with Animation */}
            {(action?.originalCode || action?.fixedCode) && (
              <div className="relative">
                <div className="flex items-center gap-2 mb-1">
                  <Code2 className="w-3 h-3 text-gray-500" />
                  <span className="text-[10px] text-gray-500 font-medium uppercase">
                    {action.action === 'healing' || action.action === 'healed'
                      ? 'Fixed Code'
                      : 'Original Code'}
                  </span>
                  {isTyping && (
                    <Badge variant="secondary" className="text-[9px] animate-pulse">
                      Typing...
                    </Badge>
                  )}
                </div>

                <ScrollArea className="h-[150px] rounded-md border bg-slate-900">
                  <pre
                    ref={codeRef}
                    className="p-3 text-[11px] font-mono text-green-400 whitespace-pre-wrap"
                  >
                    {action.action === 'healing' || action.action === 'healed'
                      ? displayedCode || action.fixedCode
                      : action.originalCode}
                    {isTyping && (
                      <motion.span
                        animate={{ opacity: [0, 1] }}
                        transition={{ repeat: Infinity, duration: 0.5 }}
                        className="inline-block w-2 h-4 bg-green-400 ml-0.5"
                      />
                    )}
                  </pre>

                  {/* Error underline animation for original code */}
                  {action.action === 'detecting' && (
                    <motion.div
                      initial={{ width: 0 }}
                      animate={{ width: '60%' }}
                      className="absolute bottom-8 left-3 h-0.5 bg-red-500"
                    />
                  )}
                </ScrollArea>
              </div>
            )}

            {/* Fix Explanation */}
            {action?.fixExplanation && (action.action === 'healing' || action.action === 'healed') && (
              <motion.div
                initial={{ opacity: 0, y: 10 }}
                animate={{ opacity: 1, y: 0 }}
                className="p-2 bg-green-100 border border-green-300 rounded-md"
              >
                <p className="text-[10px] text-green-600 font-medium mb-1">
                  <Wrench className="w-3 h-3 inline mr-1" />
                  Fix Applied
                </p>
                <p className="text-[11px] text-green-800">
                  {action.fixExplanation}
                </p>
              </motion.div>
            )}

            {/* Success Message */}
            {action?.action === 'healed' && (
              <motion.div
                initial={{ scale: 0.8, opacity: 0 }}
                animate={{ scale: 1, opacity: 1 }}
                className="flex items-center justify-center gap-2 p-3 bg-green-100 rounded-lg"
              >
                <CheckCircle2 className="w-5 h-5 text-green-600" />
                <span className="text-sm font-medium text-green-700">
                  Code Healed Successfully!
                </span>
              </motion.div>
            )}
          </CardContent>
        </Card>
      </motion.div>
    </AnimatePresence>
  )
}

/**
 * Compact healer status indicator
 */
interface HealerStatusIndicatorProps {
  isActive: boolean
  status: HealerAction['action'] | null
}

export function HealerStatusIndicator({
  isActive,
  status
}: HealerStatusIndicatorProps) {
  if (!isActive) return null

  const getConfig = () => {
    switch (status) {
      case 'detecting':
        return { color: 'bg-red-500', label: 'Error Detected' }
      case 'analyzing':
        return { color: 'bg-amber-500', label: 'Analyzing' }
      case 'healing':
        return { color: 'bg-blue-500', label: 'Healing' }
      case 'healed':
        return { color: 'bg-green-500', label: 'Healed' }
      case 'failed':
        return { color: 'bg-red-500', label: 'Failed' }
      default:
        return { color: 'bg-gray-500', label: 'Monitoring' }
    }
  }

  const config = getConfig()

  return (
    <Badge variant="outline" className="gap-1 text-xs">
      <motion.span
        animate={{ scale: [1, 1.3, 1] }}
        transition={{ repeat: Infinity, duration: 1 }}
        className={`w-2 h-2 rounded-full ${config.color}`}
      />
      {config.label}
    </Badge>
  )
}
