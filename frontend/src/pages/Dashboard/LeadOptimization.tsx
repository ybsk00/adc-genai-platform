/**
 * Lead Optimization Page
 * 기존 ADC 후보 최적화 페이지 - Trinity Layout 적용
 */
import { useState, useCallback } from 'react'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Label } from '@/components/ui/label'
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select'
import { Textarea } from '@/components/ui/textarea'
import { Input } from '@/components/ui/input'
import { Badge } from '@/components/ui/badge'
import { Separator } from '@/components/ui/separator'
import { toast } from 'sonner'
import { motion, AnimatePresence } from 'framer-motion'
import {
  Sparkles,
  Play,
  Loader2,
  Lock,
  CheckCircle2,
  AlertTriangle,
  FlaskConical,
  TrendingUp,
  ChevronDown,
  ChevronUp
} from 'lucide-react'
import { useDesignSession, type Candidate } from '@/hooks/useDesignSession'

// Design Components
import { TrinityTabsLayout } from '@/components/design/TrinityLayout'
import { LiveAgentConsole } from '@/components/design/LiveAgentConsole'
import { DesignProgressTimeline } from '@/components/design/DesignProgressTimeline'
import { MoleculeDiffViewer } from '@/components/design/MoleculeViewer2D'
import { DeltaMetricDisplay, createOptimizationMetrics, type MetricDelta } from '@/components/design/DeltaMetricDisplay'
import { EvidenceHub, type LiteratureReference, type GoldenSetReference } from '@/components/design/EvidenceHub'
import { HealerAnimation } from '@/components/design/HealerAnimation'
import { DigitalSealBadge } from '@/components/design/DigitalSealBadge'

const API_BASE_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000'

interface OptimizationFormData {
  source_smiles: string
  optimization_target: string
  target_property: string
  property_direction: 'increase' | 'decrease' | 'optimize'
  constraints: string
}

const OPTIMIZATION_TARGETS = [
  { value: 'potency', label: 'Potency (IC50)' },
  { value: 'selectivity', label: 'Selectivity' },
  { value: 'stability', label: 'Stability' },
  { value: 'solubility', label: 'Solubility' },
  { value: 'permeability', label: 'Permeability' },
  { value: 'metabolic_stability', label: 'Metabolic Stability' }
]

const PROPERTY_DIRECTIONS = [
  { value: 'increase', label: 'Increase' },
  { value: 'decrease', label: 'Decrease' },
  { value: 'optimize', label: 'Optimize (AI decides)' }
]

export default function LeadOptimization() {
  const [sessionId, setSessionId] = useState<string | null>(null)
  const [isCreating, setIsCreating] = useState(false)
  const [showAdvanced, setShowAdvanced] = useState(false)
  const [selectedCandidateIndex, setSelectedCandidateIndex] = useState<number>(0)
  const [formData, setFormData] = useState<OptimizationFormData>({
    source_smiles: '',
    optimization_target: 'potency',
    target_property: '',
    property_direction: 'optimize',
    constraints: ''
  })

  const {
    isConnected,
    agentLogs,
    currentAgent,
    currentStep,
    candidates,
    sessionStatus,
    // New state
    calculatedMetrics,
    healerAction,
    isHealerActive,
    literatureRefs,
    goldenSetRefs,
    digitalSeal
  } = useDesignSession({
    sessionId: sessionId || '',
    onComplete: (event) => {
      if (event.status === 'completed') {
        toast.success('Optimization completed successfully!')
      } else if (event.status === 'manual_review') {
        toast.warning('Optimization requires manual review')
      }
    },
    onError: (error) => {
      toast.error(`Connection error: ${error.message}`)
    }
  })

  const handleInputChange = (field: keyof OptimizationFormData, value: string) => {
    setFormData(prev => ({ ...prev, [field]: value }))
  }

  const handleCreateSession = useCallback(async () => {
    if (!formData.source_smiles) {
      toast.error('Please enter a source SMILES')
      return
    }

    setIsCreating(true)

    try {
      const createResponse = await fetch(`${API_BASE_URL}/api/design/session`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          session_type: 'optimization',
          design_goal: `Optimize ${formData.optimization_target}: ${formData.property_direction}. Source: ${formData.source_smiles}. ${formData.constraints}`
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
        throw new Error('Failed to start optimization')
      }

      toast.success('Optimization started! Watch the progress below.')
    } catch (error) {
      console.error('Optimization error:', error)
      toast.error('Failed to start optimization')
    } finally {
      setIsCreating(false)
    }
  }, [formData])

  const isRunning = sessionStatus === 'running' || Boolean(sessionId && !sessionStatus)

  // Get selected optimized candidate
  const selectedCandidate = candidates[selectedCandidateIndex]

  // Create optimization metrics for DeltaMetricDisplay
  const optimizationMetrics: MetricDelta[] = calculatedMetrics
    ? createOptimizationMetrics({
        toxicityBefore: calculatedMetrics.original_toxicity,
        toxicityAfter: calculatedMetrics.optimized_toxicity,
        stabilityBefore: calculatedMetrics.original_stability,
        stabilityAfter: calculatedMetrics.optimized_stability,
        saScoreBefore: calculatedMetrics.original_sa_score,
        saScoreAfter: calculatedMetrics.optimized_sa_score,
        logPBefore: calculatedMetrics.original_logP,
        logPAfter: calculatedMetrics.optimized_logP
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
      <Card>
        <CardHeader>
          <CardTitle className="text-base flex items-center gap-2">
            <TrendingUp className="w-4 h-4" />
            Optimization Parameters
          </CardTitle>
          <CardDescription>
            Define your lead compound and optimization goals
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          {/* Source SMILES */}
          <div className="space-y-2">
            <Label htmlFor="source_smiles">Source SMILES *</Label>
            <Textarea
              id="source_smiles"
              placeholder="Enter the SMILES of your lead compound..."
              value={formData.source_smiles}
              onChange={(e) => handleInputChange('source_smiles', e.target.value)}
              disabled={isRunning}
              rows={3}
              className="font-mono text-sm"
            />
          </div>

          {/* Optimization Target */}
          <div className="space-y-2">
            <Label htmlFor="optimization_target">Optimization Target</Label>
            <Select
              value={formData.optimization_target}
              onValueChange={(value) => handleInputChange('optimization_target', value)}
              disabled={isRunning}
            >
              <SelectTrigger>
                <SelectValue placeholder="Select target property" />
              </SelectTrigger>
              <SelectContent>
                {OPTIMIZATION_TARGETS.map((target) => (
                  <SelectItem key={target.value} value={target.value}>
                    {target.label}
                  </SelectItem>
                ))}
              </SelectContent>
            </Select>
          </div>

          {/* Property Direction */}
          <div className="space-y-2">
            <Label htmlFor="property_direction">Direction</Label>
            <Select
              value={formData.property_direction}
              onValueChange={(value) => handleInputChange('property_direction', value as 'increase' | 'decrease' | 'optimize')}
              disabled={isRunning}
            >
              <SelectTrigger>
                <SelectValue placeholder="Select direction" />
              </SelectTrigger>
              <SelectContent>
                {PROPERTY_DIRECTIONS.map((dir) => (
                  <SelectItem key={dir.value} value={dir.value}>
                    {dir.label}
                  </SelectItem>
                ))}
              </SelectContent>
            </Select>
          </div>

          {/* Advanced Options Toggle */}
          <button
            type="button"
            onClick={() => setShowAdvanced(!showAdvanced)}
            className="flex items-center gap-1 text-sm text-gray-500 hover:text-gray-700"
          >
            {showAdvanced ? <ChevronUp className="w-4 h-4" /> : <ChevronDown className="w-4 h-4" />}
            Advanced Options
          </button>

          {showAdvanced && (
            <>
              <div className="space-y-2">
                <Label htmlFor="target_property">Target Value (Optional)</Label>
                <Input
                  id="target_property"
                  placeholder="e.g., IC50 < 10nM"
                  value={formData.target_property}
                  onChange={(e) => handleInputChange('target_property', e.target.value)}
                  disabled={isRunning}
                />
              </div>

              <div className="space-y-2">
                <Label htmlFor="constraints">Additional Constraints</Label>
                <Textarea
                  id="constraints"
                  placeholder="e.g., Maintain MW < 1000, preserve core scaffold..."
                  value={formData.constraints}
                  onChange={(e) => handleInputChange('constraints', e.target.value)}
                  disabled={isRunning}
                  rows={2}
                />
              </div>
            </>
          )}

          <Separator />

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
                Optimization in Progress...
              </>
            ) : (
              <>
                <Play className="w-4 h-4 mr-2" />
                Start Optimization
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
  )

  // === CENTER PANEL: Visualization ===
  const centerPanel = (
    <div className="space-y-4">
      {/* Digital Seal */}
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

      {/* Molecule Diff Viewer - Original vs Optimized */}
      {formData.source_smiles && selectedCandidate && !selectedCandidate.is_masked && (
        <Card>
          <CardHeader>
            <CardTitle className="text-sm font-medium">Structure Comparison</CardTitle>
          </CardHeader>
          <CardContent>
            <MoleculeDiffViewer
              originalSmiles={formData.source_smiles}
              optimizedSmiles={selectedCandidate.smiles}
              originalLabel="Original"
              optimizedLabel="Optimized"
              toxicityHighlights={[
                { id: 'tox1', atomIndices: [], color: '#ef4444', label: 'Toxicity Site', type: 'toxicity' }
              ]}
              improvementHighlights={[
                { id: 'imp1', atomIndices: [], color: '#22c55e', label: 'Improved', type: 'modification' }
              ]}
            />
          </CardContent>
        </Card>
      )}

      {/* Delta Metrics Display */}
      {optimizationMetrics.length > 0 && (
        <DeltaMetricDisplay
          metrics={optimizationMetrics}
          title="Optimization Results"
          showSummary={true}
        />
      )}

      {/* Candidates List */}
      {candidates.length > 0 && (
        <Card>
          <CardHeader>
            <CardTitle className="text-base flex items-center gap-2">
              <FlaskConical className="w-4 h-4" />
              Optimized Candidates
              <Badge variant="secondary" className="ml-auto">
                {candidates.length} variants
              </Badge>
            </CardTitle>
          </CardHeader>
          <CardContent>
            <div className="space-y-3">
              <AnimatePresence>
                {candidates.map((candidate, index) => (
                  <OptimizedCandidateCard
                    key={`candidate-${index}`}
                    candidate={candidate}
                    index={index}
                    isSelected={index === selectedCandidateIndex}
                    onSelect={() => setSelectedCandidateIndex(index)}
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
          className="flex items-center gap-2 p-4 bg-green-50 border border-green-200 rounded-lg"
        >
          <CheckCircle2 className="w-5 h-5 text-green-500" />
          <span className="text-green-700 font-medium">
            Optimization completed! Review the improved candidates above.
          </span>
        </motion.div>
      )}

      {sessionStatus === 'manual_review' && (
        <motion.div
          initial={{ opacity: 0, y: 10 }}
          animate={{ opacity: 1, y: 0 }}
          className="flex items-center gap-2 p-4 bg-amber-50 border border-amber-200 rounded-lg"
        >
          <AlertTriangle className="w-5 h-5 text-amber-500" />
          <span className="text-amber-700 font-medium">
            Optimization results require expert validation.
          </span>
        </motion.div>
      )}
    </div>
  )

  // === RIGHT PANEL: Evidence Hub ===
  const rightPanel = (
    <div className="space-y-4">
      <LiveAgentConsole
        logs={agentLogs}
        currentAgent={currentAgent}
        isConnected={isConnected}
      />

      <EvidenceHub
        references={evidenceRefs}
        goldenSetRefs={goldenRefs}
        isLoading={currentAgent === 'librarian'}
      />
    </div>
  )

  return (
    <div className="p-6 space-y-6 h-full">
      {/* Header */}
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-2xl font-bold flex items-center gap-2">
            <Sparkles className="w-6 h-6 text-amber-500" />
            Lead Optimization
          </h1>
          <p className="text-gray-500 text-sm mt-1">
            AI-powered optimization of existing ADC candidates
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
        centerLabel="Comparison"
        rightLabel="Evidence"
      />
    </div>
  )
}

// Optimized Candidate Card Component
function OptimizedCandidateCard({
  candidate,
  index,
  isSelected,
  onSelect
}: {
  candidate: Candidate
  index: number
  isSelected: boolean
  onSelect: () => void
}) {
  return (
    <motion.div
      initial={{ opacity: 0, x: -20 }}
      animate={{ opacity: 1, x: 0 }}
      transition={{ delay: index * 0.1 }}
      onClick={onSelect}
      className={`p-4 rounded-lg border cursor-pointer transition-colors ${
        candidate.is_masked
          ? 'bg-gray-50 border-gray-200'
          : isSelected
          ? 'bg-amber-50 border-amber-300 ring-2 ring-amber-200'
          : 'bg-white border-slate-200 hover:border-amber-200'
      }`}
    >
      <div className="flex items-start justify-between">
        <div className="flex items-center gap-2">
          <Badge variant="outline" className="text-xs">
            Variant #{candidate.rank}
          </Badge>
          {candidate.is_masked && (
            <Badge variant="secondary" className="text-xs flex items-center gap-1">
              <Lock className="w-3 h-3" />
              Premium Only
            </Badge>
          )}
        </div>
        {!candidate.is_masked && candidate.score && (
          <Badge className="bg-amber-500">
            Improvement: +{typeof candidate.score === 'number' ? (candidate.score * 100).toFixed(1) : candidate.score}%
          </Badge>
        )}
      </div>

      <div className="mt-3">
        <Label className="text-xs text-gray-500">SMILES</Label>
        <p
          className={`font-mono text-sm mt-1 break-all ${
            candidate.is_masked ? 'text-gray-400' : 'text-gray-700'
          }`}
        >
          {candidate.smiles}
        </p>
      </div>

      {!candidate.is_masked && candidate.metrics && (
        <div className="grid grid-cols-4 gap-2 mt-3 pt-3 border-t">
          {Object.entries(candidate.metrics).map(([key, value]) => (
            <div key={key} className="text-center">
              <p className="text-xs text-gray-500 uppercase">{key}</p>
              <p className="text-sm font-medium">
                {typeof value === 'number' ? value.toFixed(2) : value}
              </p>
            </div>
          ))}
        </div>
      )}
    </motion.div>
  )
}
