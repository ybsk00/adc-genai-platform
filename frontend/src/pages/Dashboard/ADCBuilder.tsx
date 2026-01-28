import { useState, useEffect } from 'react'
import { motion, AnimatePresence } from 'framer-motion'
import { useNavigate } from 'react-router-dom'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Textarea } from '@/components/ui/textarea'
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select'
import { Badge } from '@/components/ui/badge'
import { Skeleton } from '@/components/ui/skeleton'
import {
    ChevronRight,
    ChevronLeft,
    Dna,
    FlaskConical,
    Settings2,
    Rocket,
    Upload,
    Check,
    AlertCircle,
    Loader2,
    Database
} from 'lucide-react'
import { toast } from 'sonner'
import { useADCBuilderStore } from '@/stores/adcBuilderStore'
import { MoleculeViewer2D } from '@/components/design/MoleculeViewer2D'

const API_BASE_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000'

const steps = [
    { id: 1, title: 'Target & Antibody', icon: Dna },
    { id: 2, title: 'Payload & Linker', icon: FlaskConical },
    { id: 3, title: 'Configuration', icon: Settings2 },
]

// Types for API data
interface AntibodyOption {
    id: string
    drug_name: string
    target_1: string
    target_2?: string
    antibody_type?: string
}

interface PayloadOption {
    id: string
    name: string
    smiles: string
    category: string
    mw?: number
}

interface LinkerOption {
    id: string
    name: string
    linker_type: string
    description?: string
}

export function ADCBuilder() {
    const navigate = useNavigate()
    const [currentStep, setCurrentStep] = useState(1)
    const [isSubmitting, setIsSubmitting] = useState(false)

    // Data loading states
    const [isLoadingAntibodies, setIsLoadingAntibodies] = useState(true)
    const [isLoadingPayloads, setIsLoadingPayloads] = useState(true)
    const [isLoadingLinkers, setIsLoadingLinkers] = useState(true)

    // API data
    const [antibodyOptions, setAntibodyOptions] = useState<AntibodyOption[]>([])
    const [payloadOptions, setPayloadOptions] = useState<PayloadOption[]>([])
    const [linkerOptions, setLinkerOptions] = useState<LinkerOption[]>([])

    // Zustand store
    const {
        antibodyType, customSequence, targetName,
        payloadId, linkerId, dar,
        mode, jobName,
        sequenceError,
        setField, validateSequence, resetForm
    } = useADCBuilderStore()

    // Fetch antibodies from Golden Set
    useEffect(() => {
        const fetchAntibodies = async () => {
            try {
                const response = await fetch(`${API_BASE_URL}/api/library/golden-set?limit=20`)
                if (response.ok) {
                    const data = await response.json()
                    setAntibodyOptions([
                        { id: 'custom', drug_name: 'Custom (Manual Input)', target_1: '' },
                        ...data.items.map((item: any) => ({
                            id: item.id,
                            drug_name: item.drug_name,
                            target_1: item.target_1,
                            target_2: item.target_2,
                            antibody_type: item.antibody_type
                        }))
                    ])
                }
            } catch (error) {
                console.error('Failed to fetch antibodies:', error)
                // Fallback to basic options
                setAntibodyOptions([
                    { id: 'custom', drug_name: 'Custom (Manual Input)', target_1: '' }
                ])
            } finally {
                setIsLoadingAntibodies(false)
            }
        }
        fetchAntibodies()
    }, [])

    // Fetch payloads from commercial_reagents
    useEffect(() => {
        const fetchPayloads = async () => {
            try {
                const response = await fetch(`${API_BASE_URL}/api/library/reagents?category=payload&limit=20`)
                if (response.ok) {
                    const data = await response.json()
                    setPayloadOptions(data.items.map((item: any) => ({
                        id: item.id,
                        name: item.name,
                        smiles: item.smiles,
                        category: item.sub_category || 'Payload',
                        mw: item.molecular_weight
                    })))
                }
            } catch (error) {
                console.error('Failed to fetch payloads:', error)
                // Fallback to common payloads
                setPayloadOptions([
                    { id: 'mmae', name: 'MMAE', smiles: 'CC(C)C[C@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)c1ccc(cc1)N(C)C)C(C)C)C(=O)N[C@@H](C(C)C)C(=O)N(C)[C@@H](Cc1c[nH]c2ccccc12)C(=O)O', category: 'Auristatin', mw: 718 },
                    { id: 'dxd', name: 'DXd (Deruxtecan)', smiles: 'COc1cc2c(cc1OC)-c1cc3c(c(=O)n1C2)COC(=O)[C@]3(CC)O', category: 'Camptothecin', mw: 494 },
                    { id: 'sn38', name: 'SN-38', smiles: 'CCc1c2c(nc3ccc(O)cc13)-c1cc3c(c(=O)n1C2)COC(=O)[C@]3(CC)O', category: 'Camptothecin', mw: 392 }
                ])
            } finally {
                setIsLoadingPayloads(false)
            }
        }
        fetchPayloads()
    }, [])

    // Fetch linkers from commercial_reagents
    useEffect(() => {
        const fetchLinkers = async () => {
            try {
                const response = await fetch(`${API_BASE_URL}/api/library/reagents?category=linker&limit=20`)
                if (response.ok) {
                    const data = await response.json()
                    setLinkerOptions(data.items.map((item: any) => ({
                        id: item.id,
                        name: item.name,
                        linker_type: item.linker_type || 'Unknown',
                        description: item.description
                    })))
                }
            } catch (error) {
                console.error('Failed to fetch linkers:', error)
                // Fallback to common linkers
                setLinkerOptions([
                    { id: 'val-cit', name: 'Val-Cit-PABC', linker_type: 'Cleavable', description: 'Protease sensitive' },
                    { id: 'mcc', name: 'MCC (SMCC)', linker_type: 'Non-cleavable', description: 'Thioether linkage' },
                    { id: 'ggfg', name: 'GGFG', linker_type: 'Cleavable', description: 'Next-gen peptide linker' }
                ])
            } finally {
                setIsLoadingLinkers(false)
            }
        }
        fetchLinkers()
    }, [])

    const handleNext = () => {
        if (currentStep === 1) {
            if (antibodyType === 'custom' && sequenceError) {
                return
            }
        }
        if (currentStep < 3) setCurrentStep(currentStep + 1)
    }

    const handleBack = () => {
        if (currentStep > 1) setCurrentStep(currentStep - 1)
    }

    const handleSubmit = async () => {
        // Validation
        if (!targetName && antibodyType === 'custom') {
            toast.error('Please enter a target name')
            return
        }

        setIsSubmitting(true)

        try {
            // 1. Create design session via API
            const selectedAntibody = antibodyOptions.find(a => a.id === antibodyType)
            const selectedPayload = payloadOptions.find(p => p.id === payloadId)
            const selectedLinker = linkerOptions.find(l => l.id === linkerId)

            const response = await fetch(`${API_BASE_URL}/api/design/session`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    session_type: 'denovo',
                    target_antigen: targetName || selectedAntibody?.target_1 || '',
                    target_indication: '',
                    requested_dar: dar,
                    linker_preference: selectedLinker?.linker_type || 'any',
                    design_goal: `ADC Builder: ${selectedAntibody?.drug_name || 'Custom'} + ${selectedPayload?.name || ''} + ${selectedLinker?.name || ''}`
                })
            })

            if (!response.ok) {
                throw new Error('Failed to create session')
            }

            const { session_id, status } = await response.json()

            // 2. Start the design workflow
            const startResponse = await fetch(
                `${API_BASE_URL}/api/design/session/${session_id}/start`,
                { method: 'POST' }
            )

            if (!startResponse.ok) {
                console.warn('Start endpoint not available, continuing anyway')
            }

            toast.success('Design session started!')
            resetForm()

            // 3. Navigate to result page with actual session ID
            navigate(`/dashboard/result/${session_id}`)

        } catch (error) {
            console.error('Submit error:', error)
            toast.error('Failed to start design. Please try again.')
        } finally {
            setIsSubmitting(false)
        }
    }

    // Get selected payload for preview
    const selectedPayload = payloadOptions.find(p => p.id === payloadId)

    return (
        <div className="max-w-4xl mx-auto">
            {/* Progress Steps */}
            <div className="mb-8">
                <div className="flex items-center justify-center">
                    {steps.map((step, index) => (
                        <div key={step.id} className="flex items-center">
                            <div className="flex flex-col items-center">
                                <div
                                    className={`w-12 h-12 rounded-full flex items-center justify-center transition-colors ${currentStep >= step.id
                                        ? 'bg-blue-600 text-white'
                                        : 'bg-slate-800 text-slate-500'
                                        }`}
                                >
                                    {currentStep > step.id ? (
                                        <Check className="w-5 h-5" />
                                    ) : (
                                        <step.icon className="w-5 h-5" />
                                    )}
                                </div>
                                <span className={`text-sm mt-2 ${currentStep >= step.id ? 'text-white font-medium' : 'text-slate-500'}`}>
                                    {step.title}
                                </span>
                            </div>
                            {index < steps.length - 1 && (
                                <div className={`w-24 h-0.5 mx-2 ${currentStep > step.id ? 'bg-blue-600' : 'bg-slate-800'}`} />
                            )}
                        </div>
                    ))}
                </div>
            </div>

            {/* Step Content */}
            <AnimatePresence mode="wait">
                <motion.div
                    key={currentStep}
                    initial={{ opacity: 0, x: 20 }}
                    animate={{ opacity: 1, x: 0 }}
                    exit={{ opacity: 0, x: -20 }}
                    transition={{ duration: 0.2 }}
                >
                    {currentStep === 1 && (
                        <Card className="bg-slate-900 border-slate-800">
                            <CardHeader>
                                <CardTitle className="text-white flex items-center gap-2">
                                    Step 1: Target & Antibody
                                    <Badge variant="outline" className="text-xs border-green-500 text-green-400">
                                        <Database className="w-3 h-3 mr-1" />
                                        Live Data
                                    </Badge>
                                </CardTitle>
                                <CardDescription className="text-slate-400">
                                    Select from FDA-approved ADCs or enter custom sequence
                                </CardDescription>
                            </CardHeader>
                            <CardContent className="space-y-6">
                                <div>
                                    <label className="block text-sm font-medium text-slate-300 mb-2">
                                        Target Name
                                    </label>
                                    <Input
                                        placeholder="e.g., HER2, TROP-2, LIV-1..."
                                        value={targetName}
                                        onChange={(e) => setField('targetName', e.target.value)}
                                        className="bg-slate-950 border-slate-800 text-white placeholder:text-slate-500"
                                    />
                                </div>

                                <div>
                                    <label className="block text-sm font-medium text-slate-300 mb-2">
                                        Antibody Selection (from Golden Set)
                                    </label>
                                    {isLoadingAntibodies ? (
                                        <Skeleton className="h-10 w-full bg-slate-800" />
                                    ) : (
                                        <Select
                                            value={antibodyType}
                                            onValueChange={(value) => {
                                                setField('antibodyType', value)
                                                // Auto-fill target name from Golden Set
                                                const selected = antibodyOptions.find(a => a.id === value)
                                                if (selected && selected.target_1) {
                                                    setField('targetName', selected.target_1)
                                                }
                                            }}
                                        >
                                            <SelectTrigger className="bg-slate-950 border-slate-800 text-white">
                                                <SelectValue placeholder="Select an antibody..." />
                                            </SelectTrigger>
                                            <SelectContent className="bg-slate-900 border-slate-800 text-white max-h-80">
                                                {antibodyOptions.map((option) => (
                                                    <SelectItem
                                                        key={option.id}
                                                        value={option.id}
                                                        className="focus:bg-slate-800 focus:text-white"
                                                    >
                                                        <div className="flex items-center gap-2">
                                                            <span>{option.drug_name}</span>
                                                            {option.target_1 && (
                                                                <Badge variant="outline" className="text-[10px] border-slate-600">
                                                                    {option.target_1}
                                                                </Badge>
                                                            )}
                                                        </div>
                                                    </SelectItem>
                                                ))}
                                            </SelectContent>
                                        </Select>
                                    )}
                                </div>

                                {antibodyType === 'custom' && (
                                    <div>
                                        <label className="block text-sm font-medium text-slate-300 mb-2">
                                            FASTA Sequence
                                        </label>
                                        <Textarea
                                            placeholder={`>Antibody_Heavy_Chain\nMVLQDRSMG...`}
                                            rows={6}
                                            value={customSequence}
                                            onChange={(e) => validateSequence(e.target.value)}
                                            className={`bg-slate-950 border-slate-800 text-white placeholder:text-slate-500 font-mono text-sm ${sequenceError ? 'border-red-500 focus:ring-red-500' : ''}`}
                                        />

                                        {sequenceError && (
                                            <motion.div
                                                initial={{ opacity: 0, y: -10 }}
                                                animate={{ opacity: 1, y: 0 }}
                                                className="mt-2 p-3 bg-red-500/10 border border-red-500/20 rounded-lg flex items-start gap-2"
                                            >
                                                <AlertCircle className="w-5 h-5 text-red-400 flex-shrink-0 mt-0.5" />
                                                <div>
                                                    <p className="text-sm font-medium text-red-400">{sequenceError}</p>
                                                    <p className="text-xs text-red-500/80 mt-1">
                                                        Valid Amino Acid Codes: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
                                                    </p>
                                                </div>
                                            </motion.div>
                                        )}

                                        {!sequenceError && customSequence && (
                                            <motion.div
                                                initial={{ opacity: 0 }}
                                                animate={{ opacity: 1 }}
                                                className="mt-2 p-2 bg-green-500/10 border border-green-500/20 rounded-lg flex items-center gap-2"
                                            >
                                                <Check className="w-4 h-4 text-green-400" />
                                                <span className="text-sm text-green-400">Valid sequence</span>
                                            </motion.div>
                                        )}

                                        <div className="mt-3">
                                            <Button variant="outline" size="sm" className="border-slate-700 text-slate-300 hover:bg-slate-800 hover:text-white">
                                                <Upload className="w-4 h-4 mr-2" />
                                                Upload .pdb / .fasta file
                                            </Button>
                                        </div>
                                    </div>
                                )}
                            </CardContent>
                        </Card>
                    )}

                    {currentStep === 2 && (
                        <Card className="bg-slate-900 border-slate-800">
                            <CardHeader>
                                <CardTitle className="text-white flex items-center gap-2">
                                    Step 2: Payload & Linker
                                    <Badge variant="outline" className="text-xs border-green-500 text-green-400">
                                        <Database className="w-3 h-3 mr-1" />
                                        Live Data
                                    </Badge>
                                </CardTitle>
                                <CardDescription className="text-slate-400">
                                    Select cytotoxic payload and linker chemistry
                                </CardDescription>
                            </CardHeader>
                            <CardContent className="space-y-6">
                                <div>
                                    <label className="block text-sm font-medium text-slate-300 mb-2">
                                        Payload Selection
                                    </label>
                                    {isLoadingPayloads ? (
                                        <div className="grid grid-cols-3 gap-4">
                                            {[1, 2, 3].map(i => (
                                                <Skeleton key={i} className="h-32 bg-slate-800" />
                                            ))}
                                        </div>
                                    ) : (
                                        <div className="grid grid-cols-3 gap-4">
                                            {payloadOptions.slice(0, 6).map((option) => (
                                                <div
                                                    key={option.id}
                                                    onClick={() => setField('payloadId', option.id)}
                                                    className={`p-4 rounded-xl border-2 cursor-pointer transition-all ${payloadId === option.id
                                                        ? 'border-blue-500 bg-blue-500/10 shadow-md shadow-blue-500/10'
                                                        : 'border-slate-800 hover:border-slate-700 hover:bg-slate-800/50'
                                                        }`}
                                                >
                                                    {/* 2D Structure Preview */}
                                                    {option.smiles && (
                                                        <div className="mb-2 h-16 flex items-center justify-center bg-slate-950 rounded-lg border border-slate-800 overflow-hidden">
                                                            <MoleculeViewer2D
                                                                smiles={option.smiles}
                                                                width={100}
                                                                height={60}
                                                                showControls={false}
                                                                showSmiles={false}
                                                                className="border-0"
                                                            />
                                                        </div>
                                                    )}
                                                    <p className="font-semibold text-white text-sm">{option.name}</p>
                                                    <p className="text-xs text-slate-500">{option.category}</p>
                                                    {option.mw && (
                                                        <span className="text-xs text-slate-600">MW: {option.mw}</span>
                                                    )}
                                                </div>
                                            ))}
                                        </div>
                                    )}
                                </div>

                                <div>
                                    <label className="block text-sm font-medium text-slate-300 mb-2">
                                        Linker Selection
                                    </label>
                                    {isLoadingLinkers ? (
                                        <div className="grid grid-cols-3 gap-3">
                                            {[1, 2, 3].map(i => (
                                                <Skeleton key={i} className="h-16 bg-slate-800" />
                                            ))}
                                        </div>
                                    ) : (
                                        <div className="grid grid-cols-3 gap-3">
                                            {linkerOptions.slice(0, 6).map((option) => (
                                                <div
                                                    key={option.id}
                                                    onClick={() => setField('linkerId', option.id)}
                                                    className={`p-3 rounded-lg border-2 cursor-pointer transition-all ${linkerId === option.id
                                                        ? 'border-blue-500 bg-blue-500/10'
                                                        : 'border-slate-800 hover:border-slate-700 hover:bg-slate-800/50'
                                                        }`}
                                                >
                                                    <p className="font-medium text-white text-sm">{option.name}</p>
                                                    <div className="flex items-center gap-2 mt-1">
                                                        <Badge variant="outline" className="text-[10px] border-slate-700 text-slate-400">
                                                            {option.linker_type}
                                                        </Badge>
                                                    </div>
                                                    {option.description && (
                                                        <p className="text-xs text-slate-500 mt-1">{option.description}</p>
                                                    )}
                                                </div>
                                            ))}
                                        </div>
                                    )}
                                </div>

                                <div>
                                    <label className="block text-sm font-medium text-slate-300 mb-2">
                                        DAR (Drug-to-Antibody Ratio): <span className="font-bold text-blue-400">{dar}</span>
                                    </label>
                                    <input
                                        type="range"
                                        min={1}
                                        max={8}
                                        value={dar}
                                        onChange={(e) => setField('dar', parseInt(e.target.value))}
                                        className="w-full h-2 bg-slate-800 rounded-lg appearance-none cursor-pointer accent-blue-500"
                                    />
                                    <div className="flex justify-between text-xs text-slate-500 mt-1">
                                        <span>1</span>
                                        <span className="text-blue-400 font-medium">4 (Recommended)</span>
                                        <span>8</span>
                                    </div>
                                </div>
                            </CardContent>
                        </Card>
                    )}

                    {currentStep === 3 && (
                        <Card className="bg-slate-900 border-slate-800">
                            <CardHeader>
                                <CardTitle className="text-white">Step 3: Configuration</CardTitle>
                                <CardDescription className="text-slate-400">Review and start AI design workflow</CardDescription>
                            </CardHeader>
                            <CardContent className="space-y-6">
                                <div>
                                    <label className="block text-sm font-medium text-slate-300 mb-2">
                                        Analysis Mode
                                    </label>
                                    <div className="grid grid-cols-2 gap-4">
                                        <div
                                            onClick={() => setField('mode', 'fast')}
                                            className={`p-4 rounded-lg border-2 cursor-pointer transition-all ${mode === 'fast'
                                                ? 'border-blue-500 bg-blue-500/10'
                                                : 'border-slate-800 hover:border-slate-700 hover:bg-slate-800/50'
                                                }`}
                                        >
                                            <div className="flex items-center justify-between mb-2">
                                                <span className="font-medium text-white">Fast Scan</span>
                                                <Badge variant="outline" className="border-slate-700 text-slate-400">~30s</Badge>
                                            </div>
                                            <p className="text-sm text-slate-500">Quick property calculation</p>
                                        </div>
                                        <div
                                            onClick={() => setField('mode', 'deep')}
                                            className={`p-4 rounded-lg border-2 cursor-pointer transition-all ${mode === 'deep'
                                                ? 'border-blue-500 bg-blue-500/10'
                                                : 'border-slate-800 hover:border-slate-700 hover:bg-slate-800/50'
                                                }`}
                                        >
                                            <div className="flex items-center justify-between mb-2">
                                                <span className="font-medium text-white">Deep Analysis</span>
                                                <Badge className="bg-blue-600 hover:bg-blue-700">~2min</Badge>
                                            </div>
                                            <p className="text-sm text-slate-500">Full Multi-Agent workflow</p>
                                        </div>
                                    </div>
                                </div>

                                <div>
                                    <label className="block text-sm font-medium text-slate-300 mb-2">
                                        Job Name
                                    </label>
                                    <Input
                                        placeholder="e.g., HER2_MMAE_Optimization_01"
                                        value={jobName}
                                        onChange={(e) => setField('jobName', e.target.value)}
                                        className="bg-slate-950 border-slate-800 text-white placeholder:text-slate-500"
                                    />
                                </div>

                                {/* Summary */}
                                <div className="bg-slate-950 rounded-lg p-4 border border-slate-800">
                                    <h4 className="font-medium text-white mb-3">Summary</h4>
                                    <div className="grid grid-cols-2 gap-2 text-sm">
                                        <span className="text-slate-500">Target:</span>
                                        <span className="text-slate-300">{targetName || '-'}</span>
                                        <span className="text-slate-500">Antibody:</span>
                                        <span className="text-slate-300">
                                            {antibodyOptions.find(a => a.id === antibodyType)?.drug_name || '-'}
                                        </span>
                                        <span className="text-slate-500">Payload:</span>
                                        <span className="text-slate-300">
                                            {payloadOptions.find(p => p.id === payloadId)?.name || '-'}
                                        </span>
                                        <span className="text-slate-500">Linker:</span>
                                        <span className="text-slate-300">
                                            {linkerOptions.find(l => l.id === linkerId)?.name || '-'}
                                        </span>
                                        <span className="text-slate-500">DAR:</span>
                                        <span className="text-slate-300">{dar}</span>
                                        <span className="text-slate-500">Mode:</span>
                                        <span className="text-slate-300 font-medium">
                                            {mode === 'deep' ? 'Deep Analysis (Multi-Agent)' : 'Fast Scan'}
                                        </span>
                                    </div>
                                </div>
                            </CardContent>
                        </Card>
                    )}
                </motion.div>
            </AnimatePresence>

            {/* Navigation Buttons */}
            <div className="flex justify-between mt-6">
                <Button
                    variant="outline"
                    onClick={handleBack}
                    disabled={currentStep === 1}
                    className="border-slate-700 text-slate-300 hover:bg-slate-800 hover:text-white disabled:opacity-50"
                >
                    <ChevronLeft className="w-4 h-4 mr-2" />
                    Back
                </Button>

                {currentStep < 3 ? (
                    <Button
                        onClick={handleNext}
                        className="bg-blue-600 hover:bg-blue-700 text-white"
                        disabled={currentStep === 1 && antibodyType === 'custom' && !!sequenceError}
                    >
                        Next
                        <ChevronRight className="w-4 h-4 ml-2" />
                    </Button>
                ) : (
                    <Button
                        onClick={handleSubmit}
                        disabled={isSubmitting}
                        className="bg-blue-600 hover:bg-blue-700 text-white"
                    >
                        {isSubmitting ? (
                            <>
                                <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                                Starting...
                            </>
                        ) : (
                            <>
                                <Rocket className="w-4 h-4 mr-2" />
                                Start Design
                            </>
                        )}
                    </Button>
                )}
            </div>
        </div>
    )
}
