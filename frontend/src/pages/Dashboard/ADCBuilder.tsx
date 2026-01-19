import { useState } from 'react'
import { motion, AnimatePresence } from 'framer-motion'
import { useNavigate } from 'react-router-dom'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Textarea } from '@/components/ui/textarea'
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select'
import { Badge } from '@/components/ui/badge'
import {
    ChevronRight,
    ChevronLeft,
    Dna,
    FlaskConical,
    Settings2,
    Rocket,
    Upload,
    Check,
    AlertCircle
} from 'lucide-react'
import { useADCBuilderStore } from '@/stores/adcBuilderStore'

const steps = [
    { id: 1, title: 'Target & Antibody', icon: Dna },
    { id: 2, title: 'Payload & Linker', icon: FlaskConical },
    { id: 3, title: 'Configuration', icon: Settings2 },
]

const antibodyOptions = [
    { value: 'trastuzumab', label: 'Trastuzumab (HER2)' },
    { value: 'sacituzumab', label: 'Sacituzumab (TROP-2)' },
    { value: 'enfortumab', label: 'Enfortumab (Nectin-4)' },
    { value: 'custom', label: 'Custom (Manual Input)' },
]

// Payload 옵션 - 2D 화학 구조 SVG 포함
const payloadOptions = [
    {
        value: 'mmae',
        label: 'MMAE',
        fullName: 'Monomethyl Auristatin E',
        category: 'Microtubule Inhibitors',
        mw: '717.96',
        // 간단한 2D 구조 SVG (실제로는 RDKit이나 CDK로 생성)
        structure: `<svg viewBox="0 0 100 60" class="w-full h-12">
      <circle cx="20" cy="30" r="8" fill="#3B82F6" stroke="#1E40AF" stroke-width="2"/>
      <line x1="28" y1="30" x2="42" y2="30" stroke="#94a3b8" stroke-width="2"/>
      <circle cx="50" cy="30" r="8" fill="#10B981" stroke="#047857" stroke-width="2"/>
      <line x1="58" y1="30" x2="72" y2="30" stroke="#94a3b8" stroke-width="2"/>
      <circle cx="80" cy="30" r="8" fill="#F59E0B" stroke="#B45309" stroke-width="2"/>
      <text x="50" y="55" text-anchor="middle" font-size="8" fill="#94a3b8">Auristatin</text>
    </svg>`
    },
    {
        value: 'dxd',
        label: 'DXd',
        fullName: 'Deruxtecan',
        category: 'Topo1 Inhibitors',
        mw: '493.55',
        structure: `<svg viewBox="0 0 100 60" class="w-full h-12">
      <polygon points="30,15 50,5 70,15 70,35 50,45 30,35" fill="#8B5CF6" stroke="#6D28D9" stroke-width="2"/>
      <circle cx="50" cy="25" r="6" fill="#F59E0B" stroke="#B45309" stroke-width="2"/>
      <text x="50" y="55" text-anchor="middle" font-size="8" fill="#94a3b8">Camptothecin</text>
    </svg>`
    },
    {
        value: 'sn38',
        label: 'SN-38',
        fullName: '7-Ethyl-10-hydroxycamptothecin',
        category: 'DNA Damagers',
        mw: '392.40',
        structure: `<svg viewBox="0 0 100 60" class="w-full h-12">
      <rect x="25" y="10" width="50" height="30" rx="5" fill="#EF4444" stroke="#B91C1C" stroke-width="2"/>
      <circle cx="50" cy="25" r="8" fill="#FCD34D" stroke="#F59E0B" stroke-width="2"/>
      <text x="50" y="55" text-anchor="middle" font-size="8" fill="#94a3b8">Irinotecan</text>
    </svg>`
    },
]

const linkerOptions = [
    { value: 'val-cit', label: 'Val-Cit (Cleavable)', desc: 'Protease sensitive' },
    { value: 'mcc', label: 'MCC (Non-cleavable)', desc: 'Excellent stability' },
    { value: 'ggfg', label: 'GGFG (Peptide)', desc: 'Next-gen linker' },
]

export function ADCBuilder() {
    const navigate = useNavigate()
    const [currentStep, setCurrentStep] = useState(1)

    // Zustand 스토어 사용 - Step 간 이동해도 상태 유지
    const {
        antibodyType, customSequence, targetName,
        payloadId, linkerId, dar,
        mode, jobName,
        sequenceError,
        setField, validateSequence, resetForm
    } = useADCBuilderStore()

    const handleNext = () => {
        // Step 1 유효성 검사
        if (currentStep === 1) {
            if (antibodyType === 'custom' && sequenceError) {
                return // 서열 오류 시 진행 불가
            }
        }
        if (currentStep < 3) setCurrentStep(currentStep + 1)
    }

    const handleBack = () => {
        if (currentStep > 1) setCurrentStep(currentStep - 1)
    }

    const handleSubmit = async () => {
        // TODO: API call to start simulation
        console.log('Starting simulation with:', {
            antibodyType, customSequence, targetName,
            payloadId, linkerId, dar, mode, jobName
        })
        resetForm() // 폼 초기화
        navigate('/dashboard/result/job_new')
    }

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
                                <CardTitle className="text-white">Step 1: Target & Antibody</CardTitle>
                                <CardDescription className="text-slate-400">Select an antibody for analysis or enter manually</CardDescription>
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
                                        Antibody Selection
                                    </label>
                                    <Select
                                        value={antibodyType}
                                        onValueChange={(value) => setField('antibodyType', value)}
                                    >
                                        <SelectTrigger className="bg-slate-950 border-slate-800 text-white">
                                            <SelectValue placeholder="Select an antibody..." />
                                        </SelectTrigger>
                                        <SelectContent className="bg-slate-900 border-slate-800 text-white">
                                            {antibodyOptions.map((option) => (
                                                <SelectItem key={option.value} value={option.value} className="focus:bg-slate-800 focus:text-white">
                                                    {option.label}
                                                </SelectItem>
                                            ))}
                                        </SelectContent>
                                    </Select>
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
                                            className={`bg-slate-950 border-slate-800 text-white placeholder:text-slate-500 ${sequenceError ? 'border-red-500 focus:ring-red-500' : ''}`}
                                        />

                                        {/* Real-time validation error display */}
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
                                <CardTitle className="text-white">Step 2: Payload & Linker</CardTitle>
                                <CardDescription className="text-slate-400">Select Payload and Linker</CardDescription>
                            </CardHeader>
                            <CardContent className="space-y-6">
                                <div>
                                    <label className="block text-sm font-medium text-slate-300 mb-2">
                                        Payload Selection
                                    </label>
                                    {/* Card with 2D chemical structure preview */}
                                    <div className="grid grid-cols-3 gap-4">
                                        {payloadOptions.map((option) => (
                                            <div
                                                key={option.value}
                                                onClick={() => setField('payloadId', option.value)}
                                                className={`p-4 rounded-xl border-2 cursor-pointer transition-all ${payloadId === option.value
                                                    ? 'border-blue-500 bg-blue-500/10 shadow-md shadow-blue-500/10'
                                                    : 'border-slate-800 hover:border-slate-700 hover:bg-slate-800/50'
                                                    }`}
                                            >
                                                {/* 2D Chemical Structure SVG */}
                                                <div
                                                    className="mb-3 p-2 bg-slate-950 rounded-lg border border-slate-800"
                                                    dangerouslySetInnerHTML={{ __html: option.structure }}
                                                />
                                                <p className="font-semibold text-white">{option.label}</p>
                                                <p className="text-xs text-slate-500">{option.fullName}</p>
                                                <div className="mt-2 flex items-center justify-between">
                                                    <Badge variant="outline" className="text-xs border-slate-700 text-slate-400">
                                                        {option.category.split(' ')[0]}
                                                    </Badge>
                                                    <span className="text-xs text-slate-600">MW: {option.mw}</span>
                                                </div>
                                            </div>
                                        ))}
                                    </div>
                                </div>

                                <div>
                                    <label className="block text-sm font-medium text-slate-300 mb-2">
                                        Linker Selection
                                    </label>
                                    <div className="grid grid-cols-3 gap-3">
                                        {linkerOptions.map((option) => (
                                            <div
                                                key={option.value}
                                                onClick={() => setField('linkerId', option.value)}
                                                className={`p-3 rounded-lg border-2 cursor-pointer transition-all ${linkerId === option.value
                                                    ? 'border-blue-500 bg-blue-500/10'
                                                    : 'border-slate-800 hover:border-slate-700 hover:bg-slate-800/50'
                                                    }`}
                                            >
                                                <p className="font-medium text-white text-sm">{option.label}</p>
                                                <p className="text-xs text-slate-500">{option.desc}</p>
                                            </div>
                                        ))}
                                    </div>
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
                                <CardDescription className="text-slate-400">Review simulation settings and run</CardDescription>
                            </CardHeader>
                            <CardContent className="space-y-6">
                                <div>
                                    <label className="block text-sm font-medium text-slate-300 mb-2">
                                        Simulation Mode
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
                                                <Badge variant="outline" className="border-slate-700 text-slate-400">1 Credit</Badge>
                                            </div>
                                            <p className="text-sm text-slate-500">Quickly check 3D structure only</p>
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
                                                <Badge className="bg-blue-600 hover:bg-blue-700">10 Credits</Badge>
                                            </div>
                                            <p className="text-sm text-slate-500">Includes toxicity, patent, and competitor analysis</p>
                                        </div>
                                    </div>
                                </div>

                                <div>
                                    <label className="block text-sm font-medium text-slate-300 mb-2">
                                        Job Name
                                    </label>
                                    <Input
                                        placeholder="e.g., LIV-1_MMAE_Test_01"
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
                                        <span className="text-slate-300">{antibodyType || '-'}</span>
                                        <span className="text-slate-500">Payload:</span>
                                        <span className="text-slate-300">{payloadId || '-'}</span>
                                        <span className="text-slate-500">Linker:</span>
                                        <span className="text-slate-300">{linkerId || '-'}</span>
                                        <span className="text-slate-500">DAR:</span>
                                        <span className="text-slate-300">{dar}</span>
                                        <span className="text-slate-500">Mode:</span>
                                        <span className="text-slate-300 font-medium">
                                            {mode === 'deep' ? 'Deep Analysis (10 Credits)' : 'Fast Scan (1 Credit)'}
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
                    <Button onClick={handleSubmit} className="bg-blue-600 hover:bg-blue-700 text-white">
                        <Rocket className="w-4 h-4 mr-2" />
                        Run Simulation
                    </Button>
                )}
            </div>
        </div>
    )
}
