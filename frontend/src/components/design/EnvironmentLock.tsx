/**
 * EnvironmentLock Component
 * 21 CFR Part 11 준수를 위한 환경 버전 표시 컴포넌트
 */
import { useState, useEffect } from 'react'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import { Button } from '@/components/ui/button'
import { Separator } from '@/components/ui/separator'
import {
  Dialog,
  DialogContent,
  DialogDescription,
  DialogHeader,
  DialogTitle,
  DialogTrigger
} from '@/components/ui/dialog'
import { ScrollArea } from '@/components/ui/scroll-area'
import {
  Lock,
  CheckCircle2,
  AlertCircle,
  Copy,
  Download,
  Server,
  Package,
  Cpu,
  Brain,
  Database,
  Globe
} from 'lucide-react'
import { toast } from 'sonner'
import { motion } from 'framer-motion'
import { cn } from '@/lib/utils'

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'https://adc-backend-962229188169.asia-northeast3.run.app'

interface LibraryVersion {
  name: string
  version: string
  category: string
}

interface EnvironmentInfo {
  python_version: string
  platform: string
  libraries: LibraryVersion[]
  environment: string
}

// 카테고리별 아이콘 매핑
const CATEGORY_ICONS: Record<string, React.ReactNode> = {
  core: <Server className="w-4 h-4" />,
  ai: <Brain className="w-4 h-4" />,
  chemistry: <Package className="w-4 h-4" />,
  data: <Database className="w-4 h-4" />,
  web: <Globe className="w-4 h-4" />
}

const CATEGORY_LABELS: Record<string, string> = {
  core: 'Core Framework',
  ai: 'AI & ML',
  chemistry: 'Chemistry',
  data: 'Data & Storage',
  web: 'Web & Network'
}

interface EnvironmentLockProps {
  compact?: boolean
  className?: string
}

export function EnvironmentLock({ compact = false, className }: EnvironmentLockProps) {
  const [envInfo, setEnvInfo] = useState<EnvironmentInfo | null>(null)
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    const fetchEnvironment = async () => {
      try {
        const response = await fetch(`${API_BASE_URL}/api/system/environment`)
        if (!response.ok) throw new Error('Failed to fetch environment info')
        const data = await response.json()
        setEnvInfo(data)
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Unknown error')
      } finally {
        setIsLoading(false)
      }
    }

    fetchEnvironment()
  }, [])

  const handleCopyLock = async () => {
    try {
      const response = await fetch(`${API_BASE_URL}/api/system/environment/lock`)
      if (!response.ok) throw new Error('Failed to fetch lock')
      const data = await response.json()

      await navigator.clipboard.writeText(data.lock_content)
      toast.success('Environment lock copied to clipboard')
    } catch {
      toast.error('Failed to copy environment lock')
    }
  }

  const handleDownloadLock = async () => {
    try {
      const response = await fetch(`${API_BASE_URL}/api/system/environment/lock`)
      if (!response.ok) throw new Error('Failed to fetch lock')
      const data = await response.json()

      const blob = new Blob([data.lock_content], { type: 'text/plain' })
      const url = URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      a.download = 'requirements-lock.txt'
      a.click()
      URL.revokeObjectURL(url)

      toast.success('Environment lock downloaded')
    } catch {
      toast.error('Failed to download environment lock')
    }
  }

  if (isLoading) {
    return (
      <div className={cn('flex items-center gap-2 text-sm text-muted-foreground', className)}>
        <Lock className="w-4 h-4 animate-pulse" />
        <span>Loading environment...</span>
      </div>
    )
  }

  if (error || !envInfo) {
    return (
      <div className={cn('flex items-center gap-2 text-sm text-red-500', className)}>
        <AlertCircle className="w-4 h-4" />
        <span>Environment unavailable</span>
      </div>
    )
  }

  // Compact 버전 (인라인 표시)
  if (compact) {
    return (
      <Dialog>
        <DialogTrigger asChild>
          <button
            className={cn(
              'flex items-center gap-2 text-xs px-2 py-1 rounded-md',
              'bg-gray-100 hover:bg-gray-200 transition-colors cursor-pointer',
              className
            )}
          >
            <Lock className="w-3 h-3 text-green-600" />
            <span className="font-mono">Python {envInfo.python_version}</span>
            <Badge variant="outline" className="text-[10px] px-1 py-0">
              {envInfo.environment}
            </Badge>
          </button>
        </DialogTrigger>
        <EnvironmentDialog envInfo={envInfo} onCopy={handleCopyLock} onDownload={handleDownloadLock} />
      </Dialog>
    )
  }

  // Full 버전 (카드)
  return (
    <Card className={className}>
      <CardHeader className="pb-3">
        <CardTitle className="text-sm flex items-center gap-2">
          <Lock className="w-4 h-4 text-green-600" />
          Environment Lock
          <Badge
            variant={envInfo.environment === 'production' ? 'default' : 'secondary'}
            className="ml-auto text-xs"
          >
            {envInfo.environment}
          </Badge>
        </CardTitle>
        <CardDescription className="text-xs">
          21 CFR Part 11 compliant version tracking
        </CardDescription>
      </CardHeader>
      <CardContent className="space-y-3">
        {/* Python & Platform */}
        <div className="flex items-center justify-between text-sm">
          <div className="flex items-center gap-2">
            <Cpu className="w-4 h-4 text-blue-500" />
            <span>Python</span>
          </div>
          <span className="font-mono text-muted-foreground">{envInfo.python_version}</span>
        </div>

        <div className="flex items-center justify-between text-sm">
          <div className="flex items-center gap-2">
            <Server className="w-4 h-4 text-gray-500" />
            <span>Platform</span>
          </div>
          <span className="font-mono text-muted-foreground text-xs">{envInfo.platform}</span>
        </div>

        <Separator />

        {/* Key Libraries (grouped by category) */}
        <div className="space-y-2">
          <span className="text-xs text-muted-foreground">Key Libraries</span>
          <div className="grid grid-cols-2 gap-1 text-xs">
            {envInfo.libraries.slice(0, 6).map((lib) => (
              <div key={lib.name} className="flex items-center justify-between">
                <span className="text-muted-foreground truncate">{lib.name}</span>
                <span className="font-mono">{lib.version}</span>
              </div>
            ))}
          </div>
        </div>

        <Separator />

        {/* Actions */}
        <div className="flex gap-2">
          <Button variant="outline" size="sm" className="flex-1 text-xs" onClick={handleCopyLock}>
            <Copy className="w-3 h-3 mr-1" />
            Copy
          </Button>
          <Button variant="outline" size="sm" className="flex-1 text-xs" onClick={handleDownloadLock}>
            <Download className="w-3 h-3 mr-1" />
            Export
          </Button>
        </div>
      </CardContent>
    </Card>
  )
}

// Environment Dialog (상세 보기)
function EnvironmentDialog({
  envInfo,
  onCopy,
  onDownload
}: {
  envInfo: EnvironmentInfo
  onCopy: () => void
  onDownload: () => void
}) {
  // 카테고리별 그룹화
  const groupedLibs = envInfo.libraries.reduce((acc, lib) => {
    if (!acc[lib.category]) acc[lib.category] = []
    acc[lib.category].push(lib)
    return acc
  }, {} as Record<string, LibraryVersion[]>)

  return (
    <DialogContent className="max-w-lg">
      <DialogHeader>
        <DialogTitle className="flex items-center gap-2">
          <Lock className="w-5 h-5 text-green-600" />
          Environment Lock Details
        </DialogTitle>
        <DialogDescription>
          Complete environment specification for reproducibility and compliance
        </DialogDescription>
      </DialogHeader>

      <div className="space-y-4">
        {/* System Info */}
        <div className="grid grid-cols-2 gap-4">
          <div className="p-3 rounded-lg bg-gray-50">
            <div className="flex items-center gap-2 text-sm font-medium">
              <Cpu className="w-4 h-4 text-blue-500" />
              Python Version
            </div>
            <div className="mt-1 font-mono text-lg">{envInfo.python_version}</div>
          </div>
          <div className="p-3 rounded-lg bg-gray-50">
            <div className="flex items-center gap-2 text-sm font-medium">
              <Server className="w-4 h-4 text-gray-500" />
              Platform
            </div>
            <div className="mt-1 text-sm text-muted-foreground">{envInfo.platform}</div>
          </div>
        </div>

        {/* Libraries by Category */}
        <ScrollArea className="h-[300px] rounded-lg border p-3">
          <div className="space-y-4">
            {Object.entries(groupedLibs).map(([category, libs]) => (
              <motion.div
                key={category}
                initial={{ opacity: 0, y: 10 }}
                animate={{ opacity: 1, y: 0 }}
              >
                <div className="flex items-center gap-2 mb-2">
                  {CATEGORY_ICONS[category] || <Package className="w-4 h-4" />}
                  <span className="text-sm font-medium">
                    {CATEGORY_LABELS[category] || category}
                  </span>
                  <Badge variant="secondary" className="text-xs">
                    {libs.length}
                  </Badge>
                </div>
                <div className="grid grid-cols-2 gap-2 pl-6">
                  {libs.map((lib) => (
                    <div
                      key={lib.name}
                      className="flex items-center justify-between text-sm p-2 rounded bg-gray-50"
                    >
                      <span className="truncate">{lib.name}</span>
                      <code className="text-xs bg-white px-1.5 py-0.5 rounded border">
                        {lib.version}
                      </code>
                    </div>
                  ))}
                </div>
              </motion.div>
            ))}
          </div>
        </ScrollArea>

        {/* Actions */}
        <div className="flex gap-2">
          <Button variant="outline" className="flex-1" onClick={onCopy}>
            <Copy className="w-4 h-4 mr-2" />
            Copy Lock File
          </Button>
          <Button variant="default" className="flex-1" onClick={onDownload}>
            <Download className="w-4 h-4 mr-2" />
            Download
          </Button>
        </div>

        {/* Compliance Note */}
        <div className="flex items-start gap-2 p-3 rounded-lg bg-green-50 border border-green-200 text-sm">
          <CheckCircle2 className="w-4 h-4 text-green-600 mt-0.5" />
          <div>
            <span className="font-medium text-green-800">21 CFR Part 11 Compliant</span>
            <p className="text-green-700 text-xs mt-0.5">
              Environment versions are tracked for audit trail and reproducibility requirements.
            </p>
          </div>
        </div>
      </div>
    </DialogContent>
  )
}

// Compact 인라인 배지 버전
export function EnvironmentBadge({ className }: { className?: string }) {
  const [pythonVersion, setPythonVersion] = useState<string | null>(null)

  useEffect(() => {
    fetch(`${API_BASE_URL}/api/system/environment`)
      .then(res => res.json())
      .then(data => setPythonVersion(data.python_version))
      .catch(() => setPythonVersion(null))
  }, [])

  if (!pythonVersion) return null

  return (
    <Badge variant="outline" className={cn('gap-1 font-mono text-xs', className)}>
      <Lock className="w-3 h-3 text-green-600" />
      Python {pythonVersion}
    </Badge>
  )
}

export default EnvironmentLock
