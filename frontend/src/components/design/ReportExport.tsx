/**
 * ReportExport Component
 * PDF/HTML 리포트 내보내기 컴포넌트
 */
import { useState } from 'react'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import {
  Dialog,
  DialogContent,
  DialogDescription,
  DialogHeader,
  DialogTitle,
  DialogTrigger,
  DialogFooter
} from '@/components/ui/dialog'
import { Label } from '@/components/ui/label'
import { Switch } from '@/components/ui/switch'
import { Separator } from '@/components/ui/separator'
import {
  FileText,
  Download,
  Loader2,
  FileCode,
  Eye,
  CheckCircle2,
  AlertCircle
} from 'lucide-react'
import { toast } from 'sonner'
import { cn } from '@/lib/utils'

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'https://adc-backend-962229188169.asia-northeast3.run.app'

interface ReportExportProps {
  sessionId: string
  sessionType?: string
  className?: string
  compact?: boolean
}

interface ReportPreview {
  session_id: string
  session_type: string
  status: string
  target_antigen: string
  candidates_count: number
  logs_count: number
  available_formats: string[]
}

interface ExportOptions {
  format: 'pdf' | 'html'
  includeCode: boolean
}

export function ReportExport({
  sessionId,
  sessionType = 'design',
  className,
  compact = false
}: ReportExportProps) {
  const [isOpen, setIsOpen] = useState(false)
  const [isLoading, setIsLoading] = useState(false)
  const [isExporting, setIsExporting] = useState(false)
  const [preview, setPreview] = useState<ReportPreview | null>(null)
  const [options, setOptions] = useState<ExportOptions>({
    format: 'pdf',
    includeCode: false
  })

  // 미리보기 로드
  const loadPreview = async () => {
    setIsLoading(true)
    try {
      const response = await fetch(`${API_BASE_URL}/api/report/session/${sessionId}/report/preview`)
      if (!response.ok) throw new Error('Failed to load preview')
      const data = await response.json()
      setPreview(data)
    } catch (error) {
      toast.error('Failed to load report preview')
    } finally {
      setIsLoading(false)
    }
  }

  // 리포트 내보내기
  const handleExport = async () => {
    setIsExporting(true)
    try {
      const url = `${API_BASE_URL}/api/report/session/${sessionId}/report?format=${options.format}&include_code=${options.includeCode}`
      const response = await fetch(url)

      if (!response.ok) throw new Error('Export failed')

      const contentType = response.headers.get('content-type')
      const blob = await response.blob()

      // 파일 다운로드
      const downloadUrl = URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = downloadUrl
      a.download = `report-${sessionId.slice(0, 8)}.${options.format}`
      document.body.appendChild(a)
      a.click()
      document.body.removeChild(a)
      URL.revokeObjectURL(downloadUrl)

      // PDF fallback 확인
      if (response.headers.get('X-PDF-Fallback') === 'true') {
        toast.warning('PDF generation unavailable, exported as HTML')
      } else {
        toast.success(`Report exported as ${options.format.toUpperCase()}`)
      }

      setIsOpen(false)
    } catch (error) {
      toast.error('Failed to export report')
    } finally {
      setIsExporting(false)
    }
  }

  // 새 창에서 미리보기
  const handlePreviewInNewTab = () => {
    const url = `${API_BASE_URL}/api/report/session/${sessionId}/report?format=html`
    window.open(url, '_blank')
  }

  // Compact 버전
  if (compact) {
    return (
      <Dialog open={isOpen} onOpenChange={(open) => {
        setIsOpen(open)
        if (open) loadPreview()
      }}>
        <DialogTrigger asChild>
          <Button variant="outline" size="sm" className={cn('gap-2', className)}>
            <FileText className="w-4 h-4" />
            Export Report
          </Button>
        </DialogTrigger>
        <ReportDialog
          preview={preview}
          isLoading={isLoading}
          isExporting={isExporting}
          options={options}
          setOptions={setOptions}
          onExport={handleExport}
          onPreview={handlePreviewInNewTab}
        />
      </Dialog>
    )
  }

  // Full 버전 (Card)
  return (
    <Card className={className}>
      <CardHeader className="pb-3">
        <CardTitle className="text-sm flex items-center gap-2">
          <FileText className="w-4 h-4" />
          Report Export
        </CardTitle>
        <CardDescription className="text-xs">
          Generate PDF or HTML report for this session
        </CardDescription>
      </CardHeader>
      <CardContent className="space-y-3">
        <div className="flex items-center justify-between text-sm">
          <span className="text-muted-foreground">Session Type</span>
          <Badge variant="secondary">{sessionType}</Badge>
        </div>

        <Separator />

        <div className="flex gap-2">
          <Dialog open={isOpen} onOpenChange={(open) => {
            setIsOpen(open)
            if (open) loadPreview()
          }}>
            <DialogTrigger asChild>
              <Button variant="default" size="sm" className="flex-1">
                <Download className="w-4 h-4 mr-2" />
                Export
              </Button>
            </DialogTrigger>
            <ReportDialog
              preview={preview}
              isLoading={isLoading}
              isExporting={isExporting}
              options={options}
              setOptions={setOptions}
              onExport={handleExport}
              onPreview={handlePreviewInNewTab}
            />
          </Dialog>

          <Button variant="outline" size="sm" onClick={handlePreviewInNewTab}>
            <Eye className="w-4 h-4" />
          </Button>
        </div>
      </CardContent>
    </Card>
  )
}

// 리포트 다이얼로그 컴포넌트
function ReportDialog({
  preview,
  isLoading,
  isExporting,
  options,
  setOptions,
  onExport,
  onPreview
}: {
  preview: ReportPreview | null
  isLoading: boolean
  isExporting: boolean
  options: ExportOptions
  setOptions: (options: ExportOptions) => void
  onExport: () => void
  onPreview: () => void
}) {
  return (
    <DialogContent className="sm:max-w-md">
      <DialogHeader>
        <DialogTitle className="flex items-center gap-2">
          <FileText className="w-5 h-5" />
          Export Report
        </DialogTitle>
        <DialogDescription>
          Generate a comprehensive report for this design session
        </DialogDescription>
      </DialogHeader>

      {isLoading ? (
        <div className="flex items-center justify-center py-8">
          <Loader2 className="w-8 h-8 animate-spin text-blue-500" />
        </div>
      ) : preview ? (
        <div className="space-y-4">
          {/* 세션 정보 */}
          <div className="grid grid-cols-2 gap-3 text-sm">
            <div className="p-3 bg-gray-50 rounded-lg">
              <span className="text-muted-foreground text-xs">Session Type</span>
              <div className="font-medium mt-1">{preview.session_type}</div>
            </div>
            <div className="p-3 bg-gray-50 rounded-lg">
              <span className="text-muted-foreground text-xs">Status</span>
              <div className="font-medium mt-1 flex items-center gap-1">
                {preview.status === 'completed' ? (
                  <CheckCircle2 className="w-4 h-4 text-green-500" />
                ) : (
                  <AlertCircle className="w-4 h-4 text-amber-500" />
                )}
                {preview.status}
              </div>
            </div>
            <div className="p-3 bg-gray-50 rounded-lg">
              <span className="text-muted-foreground text-xs">Candidates</span>
              <div className="font-medium mt-1">{preview.candidates_count}</div>
            </div>
            <div className="p-3 bg-gray-50 rounded-lg">
              <span className="text-muted-foreground text-xs">Log Entries</span>
              <div className="font-medium mt-1">{preview.logs_count}</div>
            </div>
          </div>

          <Separator />

          {/* 포맷 선택 */}
          <div className="space-y-3">
            <Label className="text-sm font-medium">Export Format</Label>
            <div className="flex gap-2">
              <Button
                variant={options.format === 'pdf' ? 'default' : 'outline'}
                size="sm"
                className="flex-1"
                onClick={() => setOptions({ ...options, format: 'pdf' })}
              >
                <FileText className="w-4 h-4 mr-2" />
                PDF
              </Button>
              <Button
                variant={options.format === 'html' ? 'default' : 'outline'}
                size="sm"
                className="flex-1"
                onClick={() => setOptions({ ...options, format: 'html' })}
              >
                <FileCode className="w-4 h-4 mr-2" />
                HTML
              </Button>
            </div>
          </div>

          {/* 옵션 */}
          <div className="flex items-center justify-between">
            <div className="space-y-0.5">
              <Label htmlFor="include-code" className="text-sm">Include Code</Label>
              <p className="text-xs text-muted-foreground">
                Include executed Python code in report
              </p>
            </div>
            <Switch
              id="include-code"
              checked={options.includeCode}
              onCheckedChange={(checked) => setOptions({ ...options, includeCode: checked })}
            />
          </div>
        </div>
      ) : (
        <div className="text-center py-8 text-muted-foreground">
          Failed to load preview
        </div>
      )}

      <DialogFooter className="flex gap-2 sm:gap-2">
        <Button variant="outline" onClick={onPreview} disabled={isLoading}>
          <Eye className="w-4 h-4 mr-2" />
          Preview
        </Button>
        <Button onClick={onExport} disabled={isLoading || isExporting}>
          {isExporting ? (
            <>
              <Loader2 className="w-4 h-4 mr-2 animate-spin" />
              Exporting...
            </>
          ) : (
            <>
              <Download className="w-4 h-4 mr-2" />
              Export {options.format.toUpperCase()}
            </>
          )}
        </Button>
      </DialogFooter>
    </DialogContent>
  )
}

export default ReportExport
