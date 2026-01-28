/**
 * CMC & Sourcing Page
 * ADC 제조 및 원료 소싱 분석 페이지
 * Trinity Layout 적용 버전
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
import { ScrollArea } from '@/components/ui/scroll-area'
import { toast } from 'sonner'
import { motion, AnimatePresence } from 'framer-motion'
import {
  Factory,
  Play,
  Loader2,
  Lock,
  CheckCircle2,
  AlertTriangle,
  Package,
  Truck,
  DollarSign,
  Building2,
  Clock,
  Star,
  Download
} from 'lucide-react'
import { LiveAgentConsole } from '@/components/design/LiveAgentConsole'
import { DesignProgressTimeline } from '@/components/design/DesignProgressTimeline'
import {
  TrinityTabsLayout,
  MoleculeViewer2D,
  BOMTable,
  SAScoreGauge,
  EvidenceHub,
  DigitalSealBadge,
  EnvironmentLock
} from '@/components/design'
import { useDesignSession } from '@/hooks/useDesignSession'

const API_BASE_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000'

interface CMCFormData {
  candidate_smiles: string
  candidate_name: string
  scale: string
  timeline: string
  budget_range: string
  requirements: string
}

const SCALE_OPTIONS = [
  { value: 'mg', label: 'Milligram Scale (R&D)' },
  { value: 'gram', label: 'Gram Scale (Preclinical)' },
  { value: 'kg', label: 'Kilogram Scale (Clinical)' },
  { value: 'commercial', label: 'Commercial Scale' }
]

const TIMELINE_OPTIONS = [
  { value: 'urgent', label: '< 4 weeks (Urgent)' },
  { value: 'standard', label: '4-8 weeks (Standard)' },
  { value: 'flexible', label: '8-12 weeks (Flexible)' },
  { value: 'planning', label: '> 12 weeks (Planning)' }
]

const BUDGET_OPTIONS = [
  { value: 'low', label: '< $10K' },
  { value: 'medium', label: '$10K - $50K' },
  { value: 'high', label: '$50K - $200K' },
  { value: 'enterprise', label: '> $200K' }
]

interface Supplier {
  id: string
  name: string
  location: string
  capability: string[]
  lead_time: string
  estimated_cost: string
  rating: number
  certifications: string[]
  is_masked?: boolean
}

interface SynthesisRoute {
  route_id: string
  name: string
  steps: number
  complexity: 'low' | 'medium' | 'high'
  estimated_yield: string
  key_challenges: string[]
  is_masked?: boolean
}

interface BOMItem {
  name: string
  cas?: string
  quantity: string
  unit: string
  unitCost: number
  totalCost: number
  supplier?: string
  leadTime?: string
  isMasked?: boolean
}

export default function CMCSourcing() {
  const [sessionId, setSessionId] = useState<string | null>(null)
  const [isCreating, setIsCreating] = useState(false)
  const [suppliers, setSuppliers] = useState<Supplier[]>([])
  const [synthesisRoutes, setSynthesisRoutes] = useState<SynthesisRoute[]>([])
  const [bomItems, setBomItems] = useState<BOMItem[]>([])
  const [saScore, setSaScore] = useState<number>(5)
  const [formData, setFormData] = useState<CMCFormData>({
    candidate_smiles: '',
    candidate_name: '',
    scale: 'gram',
    timeline: 'standard',
    budget_range: 'medium',
    requirements: ''
  })

  const {
    isConnected,
    agentLogs,
    currentAgent,
    currentStep,
    sessionStatus,
    currentSmiles,
    calculatedMetrics,
    digitalSeal,
    literatureRefs,
    goldenSetRefs
  } = useDesignSession({
    sessionId: sessionId || '',
    onComplete: (event) => {
      if (event.status === 'completed') {
        toast.success('CMC analysis completed!')
        if (event.final_report?.suppliers) {
          setSuppliers(event.final_report.suppliers as Supplier[])
        }
        if (event.final_report?.synthesis_routes) {
          setSynthesisRoutes(event.final_report.synthesis_routes as SynthesisRoute[])
        }
        if (event.final_report?.bom_items) {
          setBomItems(event.final_report.bom_items as BOMItem[])
        }
        if (event.final_report?.sa_score) {
          setSaScore(event.final_report.sa_score as number)
        }
      } else if (event.status === 'manual_review') {
        toast.warning('CMC analysis requires expert review')
      }
    },
    onError: (error) => {
      toast.error(`Connection error: ${error.message}`)
    }
  })

  const handleInputChange = (field: keyof CMCFormData, value: string) => {
    setFormData(prev => ({ ...prev, [field]: value }))
  }

  const handleCreateSession = useCallback(async () => {
    if (!formData.candidate_smiles) {
      toast.error('Please enter a candidate SMILES')
      return
    }

    setIsCreating(true)
    setSuppliers([])
    setSynthesisRoutes([])
    setBomItems([])

    try {
      const createResponse = await fetch(`${API_BASE_URL}/api/design/session`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          session_type: 'cmc',
          design_goal: `CMC & Sourcing for ${formData.candidate_name || 'candidate'}. Scale: ${formData.scale}, Timeline: ${formData.timeline}, Budget: ${formData.budget_range}. SMILES: ${formData.candidate_smiles}. ${formData.requirements}`
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
        throw new Error('Failed to start CMC analysis')
      }

      toast.success('CMC analysis started! Watch the progress below.')
    } catch (error) {
      console.error('CMC error:', error)
      toast.error('Failed to start CMC analysis')
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
      a.download = `cmc-bundle-${sessionId.slice(0, 8)}.json`
      a.click()
      URL.revokeObjectURL(url)

      toast.success('CMC bundle exported')
    } catch {
      toast.error('Failed to export CMC bundle')
    }
  }

  const isRunning = sessionStatus === 'running' || Boolean(sessionId && !sessionStatus)

  // Calculate total BOM cost
  const totalBOMCost = bomItems.reduce((sum, item) => sum + (item.totalCost || 0), 0)

  // === Trinity Layout Panels ===

  // Left Panel: Input Form
  const LeftPanel = (
    <ScrollArea className="h-full">
      <div className="space-y-4 p-4">
        <Card>
          <CardHeader className="pb-3">
            <CardTitle className="text-base flex items-center gap-2">
              <Package className="w-4 h-4" />
              Manufacturing Requirements
            </CardTitle>
            <CardDescription>
              Define your candidate and production needs
            </CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <div className="space-y-2">
              <Label htmlFor="candidate_name">Candidate Name</Label>
              <Input
                id="candidate_name"
                placeholder="e.g., ADC-001, Compound-X"
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

            <Separator />

            <div className="space-y-2">
              <Label htmlFor="scale">Production Scale</Label>
              <Select
                value={formData.scale}
                onValueChange={(value) => handleInputChange('scale', value)}
                disabled={isRunning}
              >
                <SelectTrigger>
                  <SelectValue placeholder="Select scale" />
                </SelectTrigger>
                <SelectContent>
                  {SCALE_OPTIONS.map((option) => (
                    <SelectItem key={option.value} value={option.value}>
                      {option.label}
                    </SelectItem>
                  ))}
                </SelectContent>
              </Select>
            </div>

            <div className="space-y-2">
              <Label htmlFor="timeline" className="flex items-center gap-1">
                <Clock className="w-3 h-3" />
                Timeline
              </Label>
              <Select
                value={formData.timeline}
                onValueChange={(value) => handleInputChange('timeline', value)}
                disabled={isRunning}
              >
                <SelectTrigger>
                  <SelectValue placeholder="Select timeline" />
                </SelectTrigger>
                <SelectContent>
                  {TIMELINE_OPTIONS.map((option) => (
                    <SelectItem key={option.value} value={option.value}>
                      {option.label}
                    </SelectItem>
                  ))}
                </SelectContent>
              </Select>
            </div>

            <div className="space-y-2">
              <Label htmlFor="budget_range" className="flex items-center gap-1">
                <DollarSign className="w-3 h-3" />
                Budget Range
              </Label>
              <Select
                value={formData.budget_range}
                onValueChange={(value) => handleInputChange('budget_range', value)}
                disabled={isRunning}
              >
                <SelectTrigger>
                  <SelectValue placeholder="Select budget" />
                </SelectTrigger>
                <SelectContent>
                  {BUDGET_OPTIONS.map((option) => (
                    <SelectItem key={option.value} value={option.value}>
                      {option.label}
                    </SelectItem>
                  ))}
                </SelectContent>
              </Select>
            </div>

            <div className="space-y-2">
              <Label htmlFor="requirements">Special Requirements</Label>
              <Textarea
                id="requirements"
                placeholder="e.g., GMP required, specific certifications, regional preferences..."
                value={formData.requirements}
                onChange={(e) => handleInputChange('requirements', e.target.value)}
                disabled={isRunning}
                rows={2}
              />
            </div>

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
                  Analysis in Progress...
                </>
              ) : (
                <>
                  <Play className="w-4 h-4 mr-2" />
                  Start CMC Analysis
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

  // Center Panel: Results
  const CenterPanel = (
    <ScrollArea className="h-full">
      <div className="space-y-4 p-4">
        {/* Molecule Viewer + SA Score */}
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <Card>
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
                width={280}
                height={200}
              />
            </CardContent>
          </Card>

          <Card>
            <CardHeader className="pb-2">
              <CardTitle className="text-sm">Synthetic Accessibility</CardTitle>
            </CardHeader>
            <CardContent className="flex justify-center">
              <SAScoreGauge score={calculatedMetrics?.sa_score || saScore} />
            </CardContent>
          </Card>
        </div>

        {/* BOM Table */}
        {bomItems.length > 0 && (
          <BOMTable
            items={bomItems.map((item, idx) => ({
              id: `bom-${idx}`,
              name: item.name,
              casNumber: item.cas,
              quantity: `${item.quantity} ${item.unit}`,
              supplier: item.supplier || 'TBD',
              unitPrice: item.unitCost,
              currency: 'USD',
              availability: 'in-stock' as const,
              leadTime: item.leadTime
            }))}
          />
        )}

        {/* Synthesis Routes */}
        {synthesisRoutes.length > 0 && (
          <Card>
            <CardHeader className="pb-2">
              <CardTitle className="text-base flex items-center gap-2">
                <Package className="w-4 h-4" />
                Synthesis Routes
                <Badge variant="secondary" className="ml-auto">
                  {synthesisRoutes.length} routes
                </Badge>
              </CardTitle>
            </CardHeader>
            <CardContent>
              <div className="space-y-3">
                <AnimatePresence>
                  {synthesisRoutes.map((route, index) => (
                    <SynthesisRouteCard key={route.route_id} route={route} index={index} />
                  ))}
                </AnimatePresence>
              </div>
            </CardContent>
          </Card>
        )}

        {/* Suppliers */}
        {suppliers.length > 0 && (
          <Card>
            <CardHeader className="pb-2">
              <CardTitle className="text-base flex items-center gap-2">
                <Building2 className="w-4 h-4" />
                Recommended Suppliers
                <Badge variant="secondary" className="ml-auto">
                  {suppliers.length} suppliers
                </Badge>
              </CardTitle>
            </CardHeader>
            <CardContent>
              <div className="space-y-3">
                <AnimatePresence>
                  {suppliers.map((supplier, index) => (
                    <SupplierCard key={supplier.id} supplier={supplier} index={index} />
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
        {sessionStatus === 'completed' && (
          <motion.div
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            className="flex items-center gap-2 p-4 bg-green-50 border border-green-200 rounded-lg"
          >
            <CheckCircle2 className="w-5 h-5 text-green-500" />
            <span className="text-green-700 font-medium">
              CMC analysis completed! Review synthesis routes and suppliers above.
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
              Complex synthesis detected. Expert review recommended.
            </span>
          </motion.div>
        )}
      </div>
    </ScrollArea>
  )

  // Right Panel: Evidence & Export
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
          <Card>
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
                Export CMC Bundle (JSON)
              </Button>
            </CardContent>
          </Card>
        )}

        {/* Cost Summary */}
        {bomItems.length > 0 && (
          <Card>
            <CardHeader className="pb-2">
              <CardTitle className="text-sm flex items-center gap-2">
                <DollarSign className="w-4 h-4" />
                Cost Summary
              </CardTitle>
            </CardHeader>
            <CardContent className="space-y-2">
              <div className="flex justify-between text-sm">
                <span className="text-muted-foreground">Raw Materials</span>
                <span className="font-medium">${totalBOMCost.toLocaleString()}</span>
              </div>
              <Separator />
              <div className="flex justify-between text-sm font-medium">
                <span>Total Estimated</span>
                <span className="text-lg">${totalBOMCost.toLocaleString()}</span>
              </div>
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
            <Factory className="w-6 h-6 text-emerald-500" />
            CMC & Sourcing
          </h1>
          <p className="text-gray-500 text-sm mt-1">
            Manufacturing feasibility and supplier identification for ADC candidates
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
          leftLabel="Requirements"
          centerLabel="Analysis"
          rightLabel="Evidence"
        />
      </div>
    </div>
  )
}

// Synthesis Route Card Component
function SynthesisRouteCard({ route, index }: { route: SynthesisRoute; index: number }) {
  const complexityConfig = {
    low: { color: 'bg-green-100 text-green-700', label: 'Low Complexity' },
    medium: { color: 'bg-amber-100 text-amber-700', label: 'Medium Complexity' },
    high: { color: 'bg-red-100 text-red-700', label: 'High Complexity' }
  }

  const config = complexityConfig[route.complexity]

  return (
    <motion.div
      initial={{ opacity: 0, x: -20 }}
      animate={{ opacity: 1, x: 0 }}
      transition={{ delay: index * 0.1 }}
      className={`p-4 rounded-lg border ${
        route.is_masked ? 'bg-gray-50 border-gray-200' : 'bg-white border-slate-200'
      }`}
    >
      <div className="flex items-start justify-between">
        <div className="flex items-center gap-2">
          <Badge variant="outline" className="text-xs">
            Route {index + 1}
          </Badge>
          <Badge className={`text-xs ${config.color}`}>
            {config.label}
          </Badge>
          {route.is_masked && (
            <Badge variant="secondary" className="text-xs flex items-center gap-1">
              <Lock className="w-3 h-3" />
              Premium Only
            </Badge>
          )}
        </div>
        {!route.is_masked && (
          <Badge variant="outline">
            {route.steps} steps | ~{route.estimated_yield} yield
          </Badge>
        )}
      </div>

      <div className="mt-3">
        <h4 className="font-medium">{route.name}</h4>
        {!route.is_masked && route.key_challenges.length > 0 && (
          <div className="mt-2">
            <p className="text-xs text-gray-500 mb-1">Key Challenges:</p>
            <ul className="text-sm text-gray-600 list-disc list-inside">
              {route.key_challenges.map((challenge, i) => (
                <li key={i}>{challenge}</li>
              ))}
            </ul>
          </div>
        )}
      </div>
    </motion.div>
  )
}

// Supplier Card Component
function SupplierCard({ supplier, index }: { supplier: Supplier; index: number }) {
  return (
    <motion.div
      initial={{ opacity: 0, x: -20 }}
      animate={{ opacity: 1, x: 0 }}
      transition={{ delay: index * 0.1 }}
      className={`p-4 rounded-lg border ${
        supplier.is_masked ? 'bg-gray-50 border-gray-200' : 'bg-white border-slate-200'
      }`}
    >
      <div className="flex items-start justify-between">
        <div>
          <div className="flex items-center gap-2">
            <Building2 className="w-4 h-4 text-gray-400" />
            <h4 className={`font-medium ${supplier.is_masked ? 'text-gray-400' : ''}`}>
              {supplier.name}
            </h4>
            {supplier.is_masked && (
              <Badge variant="secondary" className="text-xs flex items-center gap-1">
                <Lock className="w-3 h-3" />
                Premium Only
              </Badge>
            )}
          </div>
          {!supplier.is_masked && (
            <p className="text-sm text-gray-500 mt-1">{supplier.location}</p>
          )}
        </div>
        {!supplier.is_masked && (
          <div className="flex items-center gap-1">
            {[...Array(5)].map((_, i) => (
              <Star
                key={i}
                className={`w-4 h-4 ${
                  i < supplier.rating ? 'text-amber-400 fill-amber-400' : 'text-gray-200'
                }`}
              />
            ))}
          </div>
        )}
      </div>

      {!supplier.is_masked && (
        <>
          <div className="grid grid-cols-2 gap-4 mt-3">
            <div>
              <p className="text-xs text-gray-500">Lead Time</p>
              <p className="text-sm font-medium flex items-center gap-1">
                <Truck className="w-3 h-3" />
                {supplier.lead_time}
              </p>
            </div>
            <div>
              <p className="text-xs text-gray-500">Estimated Cost</p>
              <p className="text-sm font-medium flex items-center gap-1">
                <DollarSign className="w-3 h-3" />
                {supplier.estimated_cost}
              </p>
            </div>
          </div>

          <div className="mt-3 flex flex-wrap gap-1">
            {supplier.certifications.map((cert, i) => (
              <Badge key={i} variant="outline" className="text-xs">
                {cert}
              </Badge>
            ))}
          </div>
        </>
      )}
    </motion.div>
  )
}
