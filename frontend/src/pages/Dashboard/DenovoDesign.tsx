/**
 * De novo Design Page
 * ADC 후보 신규 설계 페이지 - Trinity Layout 적용
 */
import { useState, useCallback } from 'react'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Label } from '@/components/ui/label'
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select'
import { Textarea } from '@/components/ui/textarea'
import { Badge } from '@/components/ui/badge'
import { Separator } from '@/components/ui/separator'
import { toast } from 'sonner'
import { motion, AnimatePresence } from 'framer-motion'
import {
  Beaker,
  Play,
  Loader2,
  Lock,
  CheckCircle2,
  AlertTriangle,
  FlaskConical,
  Target,
  Link2,
  Atom,
  Activity,
  FileSearch,
  Cpu,
  ArrowRight
} from 'lucide-react'
import { useDesignSession, type Candidate } from '@/hooks/useDesignSession'

// Design Components
import { TrinityTabsLayout } from '@/components/design/TrinityLayout'
import { LiveAgentConsole } from '@/components/design/LiveAgentConsole'
import { DesignProgressTimeline } from '@/components/design/DesignProgressTimeline'
import { MoleculeViewer2D, type ScaffoldHighlight } from '@/components/design/MoleculeViewer2D'
import { SuccessRadarChart, createDefaultMetrics } from '@/components/design/SuccessRadarChart'
import { EvidenceHub, type LiteratureReference, type GoldenSetReference } from '@/components/design/EvidenceHub'
import { HealerAnimation } from '@/components/design/HealerAnimation'
import { DigitalSealBadge } from '@/components/design/DigitalSealBadge'
import { AntigenAutocomplete } from '@/components/design/AntigenAutocomplete'
import { EnvironmentLock } from '@/components/design/EnvironmentLock'

const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || 'https://adc-backend-962229188169.asia-northeast3.run.app'

interface DesignFormData {
  target_antigen: string
  target_indication: string
  requested_dar: number
  linker_preference: string
  design_goal: string
}

const LINKER_OPTIONS = [
  { value: 'any', label: 'Any (AI Recommended)' },
  { value: 'cleavable', label: 'Cleavable' },
  { value: 'non-cleavable', label: 'Non-cleavable' },
  { value: 'mc-vc-pabc', label: 'MC-VC-PABC' },
  { value: 'smcc', label: 'SMCC' }
]

const DAR_OPTIONS = [2, 4, 6, 8]

export default function DenovoDesign() {
  const [sessionId, setSessionId] = useState<string | null>(null)
  const [isCreating, setIsCreating] = useState(false)
  const [selectedScaffoldId, setSelectedScaffoldId] = useState<string | null>(null)
  const [formData, setFormData] = useState<DesignFormData>({
    target_antigen: '',
    target_indication: '',
    requested_dar: 4,
    linker_preference: 'any',
    design_goal: ''
  })

  // WebSocket connection for real-time updates
  const {
    isConnected,
    agentLogs,
    currentAgent,
    currentStep,
    candidates,
    sessionStatus,
    // New state from enhanced hook
    currentSmiles,
    calculatedMetrics,
    scaffolds,
    healerAction,
    isHealerActive,
    literatureRefs,
    goldenSetRefs,
    digitalSeal
  } = useDesignSession({
    sessionId: sessionId || '',
    onComplete: (event) => {
      if (event.status === 'completed') {
        toast.success('Design completed successfully!')
      } else if (event.status === 'manual_review') {
        toast.warning('Design requires manual review')
      }
    },
    onError: (error) => {
      toast.error(`Connection error: ${error.message}`)
    }
  })

  const handleInputChange = (field: keyof DesignFormData, value: string | number) => {
    setFormData(prev => ({ ...prev, [field]: value }))
  }

  const handleCreateSession = useCallback(async () => {
    if (!formData.target_antigen) {
      toast.error('Please enter a target antigen')
      return
    }

    setIsCreating(true)

    try {
      // 1. Create session
      const createResponse = await fetch(`${API_BASE_URL}/api/design/session`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          session_type: 'denovo',
          ...formData
        })
      })

      if (!createResponse.ok) {
        throw new Error('Failed to create session')
      }

      const { session_id } = await createResponse.json()
      setSessionId(session_id)

      // 2. Start the design workflow
      const startResponse = await fetch(
        `${API_BASE_URL}/api/design/session/${session_id}/start`,
        { method: 'POST' }
      )

      if (!startResponse.ok) {
        throw new Error('Failed to start design')
      }

      toast.success('Design started! Watch the progress below.')
    } catch (error) {
      console.error('Design error:', error)
      toast.error('Failed to start design')
    } finally {
      setIsCreating(false)
    }
  }, [formData])

  const isRunning = sessionStatus === 'running' || Boolean(sessionId && !sessionStatus)

  // Convert scaffolds to highlights
  const scaffoldHighlights: ScaffoldHighlight[] = scaffolds.map((s, idx) => ({
    id: s.id,
    atomIndices: [],
    color: ['#3b82f6', '#8b5cf6', '#10b981'][idx % 3],
    label: s.type,
    type: 'scaffold'
  }))

  // Handle scaffold click for Evidence Hub linking
  const handleScaffoldClick = (scaffold: ScaffoldHighlight) => {
    setSelectedScaffoldId(scaffold.id === selectedScaffoldId ? null : scaffold.id)
  }

  // Create radar metrics from calculated metrics
  const radarMetrics = calculatedMetrics
    ? createDefaultMetrics({
      mw: calculatedMetrics.mw,
      logP: calculatedMetrics.logP,
      darMatch: calculatedMetrics.dar_match ? calculatedMetrics.dar_match * 100 : undefined,
      goldenSetSim: calculatedMetrics.golden_set_similarity ? calculatedMetrics.golden_set_similarity * 100 : undefined,
      tpsa: calculatedMetrics.tpsa,
      hbd: calculatedMetrics.hbd
    })
    : []

  // Map references for EvidenceHub
  const evidenceRefs: LiteratureReference[] = literatureRefs.map(r => ({
    id: r.id,
    pmid: r.pmid,
    title: r.title,
    authors: r.authors,
    relevanceScore: r.relevanceScore,
    summary: r.summary,
    highlightedSentences: r.highlightedSentences,
    linkedScaffoldIds: r.linkedScaffoldIds
  }))

  const goldenRefs: GoldenSetReference[] = goldenSetRefs.map(g => ({
    id: g.id,
    drugName: g.drugName,
    target: g.target,
    indication: g.indication,
    clinicalStatus: g.clinicalStatus,
    similarity: g.similarity
  }))

  // === LEFT PANEL: Input Form ===
  const leftPanel = (
    <div className="space-y-4">
      <Card className="bg-[#0f172a] border-[#1e293b]">
        <CardHeader>
          <CardTitle className="text-base flex items-center gap-2">
            <Target className="w-4 h-4" />
            Design Parameters
          </CardTitle>
          <CardDescription>
            Define your ADC design constraints
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          {/* Target Antigen with UniProt Autocomplete */}
          <div className="space-y-2">
            <Label htmlFor="target_antigen">Target Antigen *</Label>
            <AntigenAutocomplete
              value={formData.target_antigen}
              onChange={(value) => handleInputChange('target_antigen', value)}
              disabled={isRunning}
            />
          </div>

          {/* Target Indication */}
          <div className="space-y-2">
            <Label htmlFor="target_indication">Target Indication</Label>
            <Input
              id="target_indication"
              placeholder="e.g., Breast Cancer, NSCLC"
              value={formData.target_indication}
              onChange={(e) => handleInputChange('target_indication', e.target.value)}
              disabled={isRunning}
            />
          </div>

          {/* DAR */}
          <div className="space-y-2">
            <Label htmlFor="requested_dar">Drug-to-Antibody Ratio (DAR)</Label>
            <Select
              value={String(formData.requested_dar)}
              onValueChange={(value) => handleInputChange('requested_dar', parseInt(value))}
              disabled={isRunning}
            >
              <SelectTrigger>
                <SelectValue placeholder="Select DAR" />
              </SelectTrigger>
              <SelectContent>
                {DAR_OPTIONS.map((dar) => (
                  <SelectItem key={dar} value={String(dar)}>
                    DAR {dar}
                  </SelectItem>
                ))}
              </SelectContent>
            </Select>
          </div>

          {/* Linker Preference */}
          <div className="space-y-2">
            <Label htmlFor="linker_preference" className="flex items-center gap-1">
              <Link2 className="w-3 h-3" />
              Linker Preference
            </Label>
            <Select
              value={formData.linker_preference}
              onValueChange={(value) => handleInputChange('linker_preference', value)}
              disabled={isRunning}
            >
              <SelectTrigger>
                <SelectValue placeholder="Select Linker" />
              </SelectTrigger>
              <SelectContent>
                {LINKER_OPTIONS.map((option) => (
                  <SelectItem key={option.value} value={option.value}>
                    {option.label}
                  </SelectItem>
                ))}
              </SelectContent>
            </Select>
          </div>

          {/* Design Goal */}
          <div className="space-y-2">
            <Label htmlFor="design_goal">Design Goal (Optional)</Label>
            <Textarea
              id="design_goal"
              placeholder="Additional design requirements or constraints..."
              value={formData.design_goal}
              onChange={(e) => handleInputChange('design_goal', e.target.value)}
              disabled={isRunning}
              rows={3}
            />
          </div>

          <Separator />

          {/* Start Button */}
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
                Design in Progress...
              </>
            ) : (
              <>
                <Play className="w-4 h-4 mr-2" />
                Start De novo Design
              </>
            )}
          </Button>
        </CardContent>
      </Card>

      {/* Progress Timeline */}
      {sessionId && (
        <DesignProgressTimeline
          currentStep={currentStep}
          currentAgent={currentAgent}
          status={sessionStatus}
        />
      )}
    </div>
  )

  // Check if any visualization content is available
  const hasVisualizationContent = digitalSeal || isHealerActive || currentSmiles || radarMetrics.length > 0 || candidates.length > 0 || sessionStatus === 'completed' || sessionStatus === 'manual_review'

  // === CENTER PANEL: Visualization ===
  const centerPanel = (
    <div className="space-y-4">
      {/* Empty State - Welcome Guide */}
      {!hasVisualizationContent && (
        <Card className="bg-gradient-to-br from-[#0f172a] to-[#1e1b4b] border-[#1e293b]">
          <CardContent className="p-8">
            <div className="text-center space-y-6">
              {/* Icon */}
              <div className="flex justify-center">
                <div className="w-20 h-20 rounded-full bg-purple-500/10 flex items-center justify-center">
                  <Atom className="w-10 h-10 text-purple-400 animate-pulse" />
                </div>
              </div>

              {/* Title */}
              <div>
                <h3 className="text-xl font-semibold text-white mb-2">
                  AI-Powered ADC Design
                </h3>
                <p className="text-gray-400 text-sm max-w-md mx-auto">
                  Multi-Agent system designs optimal ADC candidates.
                  Set parameters in the left panel and start the design.
                </p>
              </div>

              {/* Process Steps */}
              <div className="grid grid-cols-1 md:grid-cols-4 gap-4 mt-8">
                <div className="flex flex-col items-center p-4 rounded-lg bg-slate-800/50">
                  <div className="w-10 h-10 rounded-full bg-blue-500/20 flex items-center justify-center mb-2">
                    <Cpu className="w-5 h-5 text-blue-400" />
                  </div>
                  <span className="text-xs text-gray-400 text-center">Alchemist</span>
                  <span className="text-xs text-gray-500">Molecule Generation</span>
                </div>
                <div className="flex flex-col items-center p-4 rounded-lg bg-slate-800/50">
                  <div className="w-10 h-10 rounded-full bg-green-500/20 flex items-center justify-center mb-2">
                    <Activity className="w-5 h-5 text-green-400" />
                  </div>
                  <span className="text-xs text-gray-400 text-center">Auditor</span>
                  <span className="text-xs text-gray-500">Property Evaluation</span>
                </div>
                <div className="flex flex-col items-center p-4 rounded-lg bg-slate-800/50">
                  <div className="w-10 h-10 rounded-full bg-amber-500/20 flex items-center justify-center mb-2">
                    <FileSearch className="w-5 h-5 text-amber-400" />
                  </div>
                  <span className="text-xs text-gray-400 text-center">Librarian</span>
                  <span className="text-xs text-gray-500">Literature Search</span>
                </div>
                <div className="flex flex-col items-center p-4 rounded-lg bg-slate-800/50">
                  <div className="w-10 h-10 rounded-full bg-red-500/20 flex items-center justify-center mb-2">
                    <FlaskConical className="w-5 h-5 text-red-400" />
                  </div>
                  <span className="text-xs text-gray-400 text-center">Healer</span>
                  <span className="text-xs text-gray-500">Toxicity Removal</span>
                </div>
              </div>

              {/* What you'll see */}
              <div className="border-t border-slate-700 pt-6 mt-6">
                <p className="text-xs text-gray-500 mb-3">What you'll see after starting design</p>
                <div className="flex flex-wrap justify-center gap-2">
                  <Badge variant="outline" className="text-xs">
                    <Atom className="w-3 h-3 mr-1" />
                    2D Molecular Structure
                  </Badge>
                  <Badge variant="outline" className="text-xs">
                    <Activity className="w-3 h-3 mr-1" />
                    Lipinski Radar Chart
                  </Badge>
                  <Badge variant="outline" className="text-xs">
                    <FlaskConical className="w-3 h-3 mr-1" />
                    Candidate List
                  </Badge>
                </div>
              </div>
            </div>
          </CardContent>
        </Card>
      )}

      {/* Digital Seal Badge */}
      {digitalSeal && (
        <DigitalSealBadge
          sessionId={digitalSeal.sessionId}
          recordHash={digitalSeal.recordHash}
          previousHash={digitalSeal.previousHash}
          chainHash={digitalSeal.chainHash}
          timestamp={digitalSeal.timestamp}
          isVerified={digitalSeal.isVerified}
        />
      )}

      {/* Healer Animation */}
      {isHealerActive && (
        <HealerAnimation
          action={healerAction}
          isActive={isHealerActive}
        />
      )}

      {/* Molecule Viewer */}
      {currentSmiles && (
        <MoleculeViewer2D
          smiles={currentSmiles}
          title="Current Structure"
          highlights={scaffoldHighlights}
          onScaffoldClick={handleScaffoldClick}
        />
      )}

      {/* Success Radar Chart */}
      {radarMetrics.length > 0 && (
        <SuccessRadarChart
          metrics={radarMetrics}
          title="Design Metrics"
          showTarget={true}
        />
      )}

      {/* Candidates Results */}
      {candidates.length > 0 && (
        <Card className="bg-[#0f172a] border-[#1e293b]">
          <CardHeader>
            <CardTitle className="text-base flex items-center gap-2">
              <FlaskConical className="w-4 h-4" />
              Generated Candidates
              <Badge variant="secondary" className="ml-auto">
                {candidates.length} candidates
              </Badge>
            </CardTitle>
          </CardHeader>
          <CardContent>
            <div className="space-y-3">
              <AnimatePresence>
                {candidates.map((candidate, index) => (
                  <CandidateCard
                    key={`candidate-${index}`}
                    candidate={candidate}
                    index={index}
                  />
                ))}
              </AnimatePresence>
            </div>
          </CardContent>
        </Card>
      )}

      {/* Status Messages */}
      {sessionStatus === 'completed' && (
        <motion.div
          initial={{ opacity: 0, y: 10 }}
          animate={{ opacity: 1, y: 0 }}
          className="flex items-center gap-2 p-4 bg-green-900/20 border border-green-900/50 rounded-lg"
        >
          <CheckCircle2 className="w-5 h-5 text-green-500" />
          <span className="text-green-400 font-medium">
            Design completed successfully! Review the candidates above.
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
            This design requires manual expert review before proceeding.
          </span>
        </motion.div>
      )}
    </div>
  )

  // === RIGHT PANEL: Evidence Hub ===
  const rightPanel = (
    <div className="space-y-4">
      {/* Agent Console */}
      <LiveAgentConsole
        logs={agentLogs}
        currentAgent={currentAgent}
        isConnected={isConnected}
      />

      {/* Evidence Hub - Librarian References */}
      <EvidenceHub
        references={evidenceRefs}
        goldenSetRefs={goldenRefs}
        selectedScaffoldId={selectedScaffoldId}
        isLoading={currentAgent === 'librarian'}
      />

      {/* Environment Lock - 21 CFR Part 11 Compliance */}
      <EnvironmentLock />
    </div>
  )

  return (
    <div className="p-6 space-y-6 h-full">
      {/* Header */}
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-2xl font-bold flex items-center gap-2">
            <Beaker className="w-6 h-6 text-purple-500" />
            De novo Design
          </h1>
          <p className="text-gray-500 text-sm mt-1">
            AI-powered ADC candidate generation with multi-agent orchestration
          </p>
        </div>
        {sessionId && (
          <Badge variant="outline" className="text-xs">
            Session: {sessionId.slice(0, 8)}...
          </Badge>
        )}
      </div>

      {/* Trinity Layout */}
      <TrinityTabsLayout
        leftPanel={leftPanel}
        centerPanel={centerPanel}
        rightPanel={rightPanel}
        leftLabel="Parameters"
        centerLabel="Structure"
        rightLabel="Evidence"
      />
    </div>
  )
}

// Candidate Card Component
function CandidateCard({ candidate, index }: { candidate: Candidate; index: number }) {
  return (
    <motion.div
      initial={{ opacity: 0, x: -20 }}
      animate={{ opacity: 1, x: 0 }}
      transition={{ delay: index * 0.1 }}
      className={`p-4 rounded-lg border ${candidate.is_masked
        ? 'bg-slate-900/50 border-slate-800'
        : 'bg-[#0f172a] border-[#1e293b]'
        }`}
    >
      <div className="flex items-start justify-between">
        <div className="flex items-center gap-2">
          <Badge variant="outline" className="text-xs">
            Rank #{candidate.rank}
          </Badge>
          {candidate.is_masked && (
            <Badge variant="secondary" className="text-xs flex items-center gap-1">
              <Lock className="w-3 h-3" />
              Premium Only
            </Badge>
          )}
        </div>
        {!candidate.is_masked && candidate.score && (
          <Badge className="bg-blue-500">
            Score: {typeof candidate.score === 'number' ? candidate.score.toFixed(2) : candidate.score}
          </Badge>
        )}
      </div>

      <div className="mt-3">
        <Label className="text-xs text-slate-500">SMILES</Label>
        <p
          className={`font-mono text-sm mt-1 break-all ${candidate.is_masked ? 'text-slate-500' : 'text-slate-300'
            }`}
        >
          {candidate.smiles}
        </p>
      </div>

      {!candidate.is_masked && candidate.metrics && (
        <div className="grid grid-cols-4 gap-2 mt-3 pt-3 border-t border-slate-800">
          {Object.entries(candidate.metrics).map(([key, value]) => (
            <div key={key} className="text-center">
              <p className="text-xs text-slate-500 uppercase">{key}</p>
              <p className="text-sm font-medium text-slate-300">
                {typeof value === 'number' ? value.toFixed(2) : value}
              </p>
            </div>
          ))}
        </div>
      )}
    </motion.div>
  )
}
