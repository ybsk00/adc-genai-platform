/**
 * Pre-clinical Audit Page
 * ADC 후보의 전임상 감사 페이지
 * Trinity Layout 적용 버전
 */
import { useState, useCallback } from 'react'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Label } from '@/components/ui/label'
import { Textarea } from '@/components/ui/textarea'
import { Badge } from '@/components/ui/badge'
import { Separator } from '@/components/ui/separator'
import { Checkbox } from '@/components/ui/checkbox'
import { ScrollArea } from '@/components/ui/scroll-area'
import { toast } from 'sonner'
import { motion, AnimatePresence } from 'framer-motion'
import {
  Shield,
  Play,
  Loader2,
  CheckCircle2,
  AlertTriangle,
  FileCheck,
  XCircle,
  Info,
  Download
} from 'lucide-react'
import { LiveAgentConsole } from '@/components/design/LiveAgentConsole'
import { DesignProgressTimeline } from '@/components/design/DesignProgressTimeline'
import {
  TrinityTabsLayout,
  MoleculeViewer2D,
  AuditorFeedback,
  DigitalSealBadge,
  EvidenceHub,
  AntigenAutocomplete,
  EnvironmentLock
} from '@/components/design'
import { useDesignSession } from '@/hooks/useDesignSession'

const API_BASE_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000'

interface AuditFormData {
  candidate_smiles: string
  candidate_name: string
  target_antigen: string
  audit_checks: string[]
  additional_notes: string
}

const AUDIT_CHECK_OPTIONS = [
  { id: 'lipinski', label: 'Lipinski Rule of 5', description: 'Drug-likeness assessment' },
  { id: 'pains', label: 'PAINS Filter', description: 'Pan-Assay Interference Compounds' },
  { id: 'toxicity', label: 'Toxicity Prediction', description: 'AMES, hERG, Hepatotoxicity' },
  { id: 'admet', label: 'ADMET Properties', description: 'Absorption, Distribution, Metabolism, Excretion, Toxicity' },
  { id: 'synthetic', label: 'Synthetic Accessibility', description: 'Ease of synthesis score' },
  { id: 'stability', label: 'Chemical Stability', description: 'Stability predictions' },
  { id: 'patent', label: 'Patent Landscape', description: 'Freedom-to-operate analysis' },
  { id: 'clinical', label: 'Clinical Precedent', description: 'Similar compounds in clinical trials' }
]

interface AuditResult {
  check_id: string
  check_name: string
  status: 'pass' | 'warning' | 'fail' | 'info'
  score?: number
  message: string
  details?: Record<string, unknown>
}

export default function PreclinicalAudit() {
  const [sessionId, setSessionId] = useState<string | null>(null)
  const [isCreating, setIsCreating] = useState(false)
  const [auditResults, setAuditResults] = useState<AuditResult[]>([])
  const [formData, setFormData] = useState<AuditFormData>({
    candidate_smiles: '',
    candidate_name: '',
    target_antigen: '',
    audit_checks: ['lipinski', 'pains', 'toxicity', 'admet'],
    additional_notes: ''
  })

  const {
    isConnected,
    agentLogs,
    currentAgent,
    currentStep,
    sessionStatus,
    currentSmiles,
    calculatedMetrics,
    auditFeedback,
    digitalSeal,
    literatureRefs,
    goldenSetRefs
  } = useDesignSession({
    sessionId: sessionId || '',
    onComplete: (event) => {
      if (event.status === 'completed') {
        toast.success('Audit completed successfully!')
        if (event.final_report?.audit_results) {
          setAuditResults(event.final_report.audit_results as AuditResult[])
        }
      } else if (event.status === 'manual_review') {
        toast.warning('Audit flagged issues requiring expert review')
      }
    },
    onError: (error) => {
      toast.error(`Connection error: ${error.message}`)
    }
  })

  const handleInputChange = (field: keyof AuditFormData, value: string | string[]) => {
    setFormData(prev => ({ ...prev, [field]: value }))
  }

  const handleCheckToggle = (checkId: string, checked: boolean) => {
    setFormData(prev => ({
      ...prev,
      audit_checks: checked
        ? [...prev.audit_checks, checkId]
        : prev.audit_checks.filter(id => id !== checkId)
    }))
  }

  const handleCreateSession = useCallback(async () => {
    if (!formData.candidate_smiles) {
      toast.error('Please enter a candidate SMILES')
      return
    }
    if (formData.audit_checks.length === 0) {
      toast.error('Please select at least one audit check')
      return
    }

    setIsCreating(true)
    setAuditResults([])

    try {
      const createResponse = await fetch(`${API_BASE_URL}/api/design/session`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          session_type: 'audit',
          target_antigen: formData.target_antigen,
          design_goal: `Pre-clinical audit for ${formData.candidate_name || 'candidate'}. Checks: ${formData.audit_checks.join(', ')}. SMILES: ${formData.candidate_smiles}. ${formData.additional_notes}`
        })
      })

      if (!createResponse.ok) {
        throw new Error('Failed to create session')
      }

      const { session_id } = await createResponse.json()
      setSessionId(session_id)

      const startResponse = await fetch(
        `${API_BASE_URL}/api/design/session/${session_id}/start`,
        { method: 'POST' }
      )

      if (!startResponse.ok) {
        throw new Error('Failed to start audit')
      }

      toast.success('Audit started! Watch the progress below.')
    } catch (error) {
      console.error('Audit error:', error)
      toast.error('Failed to start audit')
    } finally {
      setIsCreating(false)
    }
  }, [formData])

  const handleExportAuditBundle = async () => {
    if (!sessionId) return

    try {
      const response = await fetch(`${API_BASE_URL}/api/design/session/${sessionId}/audit-bundle`)
      if (!response.ok) throw new Error('Failed to export')

      const bundle = await response.json()
      const blob = new Blob([JSON.stringify(bundle, null, 2)], { type: 'application/json' })
      const url = URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      a.download = `audit-bundle-${sessionId.slice(0, 8)}.json`
      a.click()
      URL.revokeObjectURL(url)

      toast.success('Audit bundle exported')
    } catch {
      toast.error('Failed to export audit bundle')
    }
  }

  const isRunning = sessionStatus === 'running' || Boolean(sessionId && !sessionStatus)

  const auditSummary = {
    total: auditResults.length,
    pass: auditResults.filter(r => r.status === 'pass').length,
    warning: auditResults.filter(r => r.status === 'warning').length,
    fail: auditResults.filter(r => r.status === 'fail').length
  }

  // Transform audit results to check items for AuditorFeedback
  const mapCategory = (checkId: string): 'lipinski' | 'pains' | 'stability' | 'toxicity' | 'dar' | 'other' => {
    if (checkId.includes('lipinski')) return 'lipinski'
    if (checkId.includes('pains')) return 'pains'
    if (checkId.includes('stability')) return 'stability'
    if (checkId.includes('toxicity') || checkId.includes('admet')) return 'toxicity'
    if (checkId.includes('dar')) return 'dar'
    return 'other'
  }

  const mapStatus = (status: string): 'pass' | 'fail' | 'warning' | 'pending' => {
    if (status === 'pass') return 'pass'
    if (status === 'fail') return 'fail'
    if (status === 'warning') return 'warning'
    return 'pending'
  }

  const auditCheckItems = auditResults.map(r => ({
    id: r.check_id,
    name: r.check_name,
    category: mapCategory(r.check_id),
    status: mapStatus(r.status),
    currentValue: r.score,
    reasoning: r.message
  }))

  // === Trinity Layout Panels ===

  // Left Panel: Input Form
  const LeftPanel = (
    <ScrollArea className="h-full">
      <div className="space-y-4 p-4">
        <Card className="bg-[#0f172a] border-[#1e293b]">
          <CardHeader className="pb-3">
            <CardTitle className="text-base flex items-center gap-2">
              <FileCheck className="w-4 h-4" />
              Audit Configuration
            </CardTitle>
            <CardDescription>
              Define the candidate and select audit checks
            </CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <div className="space-y-2">
              <Label htmlFor="candidate_name">Candidate Name</Label>
              <Input
                id="candidate_name"
                placeholder="e.g., ADC-001, Lead-A"
                value={formData.candidate_name}
                onChange={(e) => handleInputChange('candidate_name', e.target.value)}
                disabled={isRunning}
              />
            </div>

            <div className="space-y-2">
              <Label htmlFor="candidate_smiles">Candidate SMILES *</Label>
              <Textarea
                id="candidate_smiles"
                placeholder="Enter the SMILES structure..."
                value={formData.candidate_smiles}
                onChange={(e) => handleInputChange('candidate_smiles', e.target.value)}
                disabled={isRunning}
                rows={3}
                className="font-mono text-sm"
              />
            </div>

            <div className="space-y-2">
              <Label htmlFor="target_antigen">Target Antigen</Label>
              <AntigenAutocomplete
                value={formData.target_antigen}
                onChange={(value) => handleInputChange('target_antigen', value)}
                disabled={isRunning}
              />
            </div>

            <Separator />

            <div className="space-y-3">
              <Label>Audit Checks *</Label>
              {AUDIT_CHECK_OPTIONS.map((check) => (
                <div key={check.id} className="flex items-start space-x-3">
                  <Checkbox
                    id={check.id}
                    checked={formData.audit_checks.includes(check.id)}
                    onCheckedChange={(checked) => handleCheckToggle(check.id, checked as boolean)}
                    disabled={isRunning}
                  />
                  <div className="grid gap-0.5 leading-none">
                    <label htmlFor={check.id} className="text-sm font-medium cursor-pointer">
                      {check.label}
                    </label>
                    <p className="text-xs text-muted-foreground">{check.description}</p>
                  </div>
                </div>
              ))}
            </div>

            <Separator />

            <div className="space-y-2">
              <Label htmlFor="additional_notes">Additional Notes</Label>
              <Textarea
                id="additional_notes"
                placeholder="Any specific concerns or focus areas..."
                value={formData.additional_notes}
                onChange={(e) => handleInputChange('additional_notes', e.target.value)}
                disabled={isRunning}
                rows={2}
              />
            </div>

            <Button
              className="w-full"
              size="lg"
              onClick={handleCreateSession}
              disabled={isCreating || isRunning}
            >
              {isCreating ? (
                <>
                  <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                  Creating Session...
                </>
              ) : isRunning ? (
                <>
                  <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                  Audit in Progress...
                </>
              ) : (
                <>
                  <Play className="w-4 h-4 mr-2" />
                  Start Audit
                </>
              )}
            </Button>
          </CardContent>
        </Card>

        {sessionId && (
          <DesignProgressTimeline
            currentStep={currentStep}
            currentAgent={currentAgent}
            status={sessionStatus}
          />
        )}
      </div>
    </ScrollArea>
  )

  // Center Panel: Visualization + Audit Results
  const CenterPanel = (
    <ScrollArea className="h-full">
      <div className="space-y-4 p-4">
        {/* Molecule Viewer */}
        <Card className="bg-[#0f172a] border-[#1e293b]">
          <CardHeader className="pb-2">
            <CardTitle className="text-sm flex items-center justify-between">
              <span>Candidate Structure</span>
              {sessionId && (
                <Badge variant="outline" className="text-xs">
                  {sessionId.slice(0, 8)}...
                </Badge>
              )}
            </CardTitle>
          </CardHeader>
          <CardContent>
            <MoleculeViewer2D
              smiles={currentSmiles || formData.candidate_smiles}
              width={400}
              height={250}
            />
          </CardContent>
        </Card>

        {/* Auditor Feedback (Redesign Loop) */}
        {(auditFeedback || auditCheckItems.length > 0) && (
          <AuditorFeedback
            isInRedesignLoop={auditFeedback?.isInRedesignLoop || false}
            checkItems={auditFeedback?.checkItems || auditCheckItems}
            redesignRequest={auditFeedback?.redesignRequest}
          />
        )}

        {/* Detailed Audit Results */}
        {auditResults.length > 0 && (
          <Card className="bg-[#0f172a] border-[#1e293b]">
            <CardHeader className="pb-2">
              <CardTitle className="text-base flex items-center gap-2">
                <FileCheck className="w-4 h-4" />
                Detailed Audit Results
                <div className="ml-auto flex gap-2">
                  <Badge variant="default" className="bg-green-500">
                    {auditSummary.pass} Pass
                  </Badge>
                  <Badge variant="secondary" className="bg-amber-500 text-white">
                    {auditSummary.warning} Warning
                  </Badge>
                  <Badge variant="destructive">
                    {auditSummary.fail} Fail
                  </Badge>
                </div>
              </CardTitle>
            </CardHeader>
            <CardContent>
              <div className="space-y-3">
                <AnimatePresence>
                  {auditResults.map((result, index) => (
                    <AuditResultCard key={result.check_id} result={result} index={index} />
                  ))}
                </AnimatePresence>
              </div>
            </CardContent>
          </Card>
        )}

        {/* Agent Console */}
        <LiveAgentConsole
          logs={agentLogs}
          currentAgent={currentAgent}
          isConnected={isConnected}
        />

        {/* Status Messages */}
        {sessionStatus === 'completed' && auditSummary.fail === 0 && (
          <motion.div
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            className="flex items-center gap-2 p-4 bg-green-900/20 border border-green-900/50 rounded-lg"
          >
            <CheckCircle2 className="w-5 h-5 text-green-500" />
            <span className="text-green-400 font-medium">
              Audit completed successfully! All checks passed.
            </span>
          </motion.div>
        )}

        {sessionStatus === 'completed' && auditSummary.fail > 0 && (
          <motion.div
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            className="flex items-center gap-2 p-4 bg-red-900/20 border border-red-900/50 rounded-lg"
          >
            <XCircle className="w-5 h-5 text-red-500" />
            <span className="text-red-400 font-medium">
              Audit completed with {auditSummary.fail} failed check(s). Review required.
            </span>
          </motion.div>
        )}

        {sessionStatus === 'manual_review' && (
          <motion.div
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            className="flex items-center gap-2 p-4 bg-amber-900/20 border border-amber-900/50 rounded-lg"
          >
            <AlertTriangle className="w-5 h-5 text-amber-500" />
            <span className="text-amber-400 font-medium">
              Audit flagged concerns requiring expert review.
            </span>
          </motion.div>
        )}
      </div>
    </ScrollArea>
  )

  // Right Panel: Evidence Hub + Digital Seal
  const RightPanel = (
    <ScrollArea className="h-full">
      <div className="space-y-4 p-4">
        {/* Digital Seal Badge */}
        {sessionId && (
          <DigitalSealBadge
            sessionId={sessionId}
            recordHash={digitalSeal?.recordHash || ''}
            chainHash={digitalSeal?.chainHash}
            timestamp={digitalSeal?.timestamp || new Date().toISOString()}
            isVerified={digitalSeal?.isVerified || false}
            onExportBundle={handleExportAuditBundle}
          />
        )}

        {/* Evidence Hub */}
        <EvidenceHub
          references={literatureRefs}
          goldenSetRefs={goldenSetRefs}
          selectedScaffoldId={null}
        />

        {/* Quick Export */}
        {sessionStatus === 'completed' && sessionId && (
          <Card className="bg-[#0f172a] border-[#1e293b]">
            <CardHeader className="pb-2">
              <CardTitle className="text-sm">Export Options</CardTitle>
            </CardHeader>
            <CardContent className="space-y-2">
              <Button
                variant="outline"
                size="sm"
                className="w-full"
                onClick={handleExportAuditBundle}
              >
                <Download className="w-4 h-4 mr-2" />
                Export Audit Bundle (JSON)
              </Button>
            </CardContent>
          </Card>
        )}

        {/* Environment Lock - 21 CFR Part 11 Compliance */}
        <EnvironmentLock />
      </div>
    </ScrollArea>
  )

  return (
    <div className="h-full flex flex-col">
      {/* Header */}
      <div className="flex items-center justify-between px-6 py-4 border-b">
        <div>
          <h1 className="text-2xl font-bold flex items-center gap-2">
            <Shield className="w-6 h-6 text-blue-500" />
            Pre-clinical Audit
          </h1>
          <p className="text-gray-500 text-sm mt-1">
            Comprehensive safety and efficacy assessment for ADC candidates
          </p>
        </div>
        {sessionId && (
          <Badge variant="outline" className="text-xs">
            Session: {sessionId.slice(0, 8)}...
          </Badge>
        )}
      </div>

      {/* Trinity Layout */}
      <div className="flex-1 overflow-hidden">
        <TrinityTabsLayout
          leftPanel={LeftPanel}
          centerPanel={CenterPanel}
          rightPanel={RightPanel}
          leftLabel="Configuration"
          centerLabel="Audit"
          rightLabel="Evidence"
        />
      </div>
    </div>
  )
}

// Audit Result Card Component
function AuditResultCard({ result, index }: { result: AuditResult; index: number }) {
  const statusConfig = {
    pass: { icon: CheckCircle2, color: 'text-green-400', bg: 'bg-green-900/20 border-green-900/50' },
    warning: { icon: AlertTriangle, color: 'text-amber-400', bg: 'bg-amber-900/20 border-amber-900/50' },
    fail: { icon: XCircle, color: 'text-red-400', bg: 'bg-red-900/20 border-red-900/50' },
    info: { icon: Info, color: 'text-blue-400', bg: 'bg-blue-900/20 border-blue-900/50' }
  }

  const config = statusConfig[result.status]
  const StatusIcon = config.icon

  return (
    <motion.div
      initial={{ opacity: 0, x: -20 }}
      animate={{ opacity: 1, x: 0 }}
      transition={{ delay: index * 0.05 }}
      className={`p-4 rounded-lg border ${config.bg}`}
    >
      <div className="flex items-start gap-3">
        <StatusIcon className={`w-5 h-5 mt-0.5 ${config.color}`} />
        <div className="flex-1">
          <div className="flex items-center justify-between">
            <h4 className="font-medium">{result.check_name}</h4>
            {result.score !== undefined && (
              <Badge variant="outline" className="text-xs">
                Score: {result.score.toFixed(2)}
              </Badge>
            )}
          </div>
          <p className="text-sm text-gray-600 mt-1">{result.message}</p>
        </div>
      </div>
    </motion.div>
  )
}
