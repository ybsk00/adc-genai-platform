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
      <line x1="28" y1="30" x2="42" y2="30" stroke="#374151" stroke-width="2"/>
      <circle cx="50" cy="30" r="8" fill="#10B981" stroke="#047857" stroke-width="2"/>
      <line x1="58" y1="30" x2="72" y2="30" stroke="#374151" stroke-width="2"/>
      <circle cx="80" cy="30" r="8" fill="#F59E0B" stroke="#B45309" stroke-width="2"/>
      <text x="50" y="55" text-anchor="middle" font-size="8" fill="#6B7280">Auristatin</text>
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
      <text x="50" y="55" text-anchor="middle" font-size="8" fill="#6B7280">Camptothecin</text>
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
      <text x="50" y="55" text-anchor="middle" font-size="8" fill="#6B7280">Irinotecan</text>
    </svg>`
    },
]

const linkerOptions = [
    { value: 'val-cit', label: 'Val-Cit (Cleavable)', desc: '프로테아제 민감' },
    { value: 'mcc', label: 'MCC (Non-cleavable)', desc: '안정성 우수' },
    { value: 'ggfg', label: 'GGFG (Peptide)', desc: '차세대 링커' },
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
                                        ? 'bg-[#007AFF] text-white'
                                        : 'bg-gray-100 text-gray-400'
                                        }`}
                                >
                                    {currentStep > step.id ? (
                                        <Check className="w-5 h-5" />
                                    ) : (
                                        <step.icon className="w-5 h-5" />
                                    )}
                                </div>
                                <span className={`text-sm mt-2 ${currentStep >= step.id ? 'text-gray-900 font-medium' : 'text-gray-400'}`}>
                                    {step.title}
                                </span>
                            </div>
                            {index < steps.length - 1 && (
                                <div className={`w-24 h-0.5 mx-2 ${currentStep > step.id ? 'bg-[#007AFF]' : 'bg-gray-200'}`} />
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
                        <Card>
                            <CardHeader>
                                <CardTitle>Step 1: Target & Antibody</CardTitle>
                                <CardDescription>분석할 항체를 선택하거나 직접 입력하세요</CardDescription>
                            </CardHeader>
                            <CardContent className="space-y-6">
                                <div>
                                    <label className="block text-sm font-medium text-gray-700 mb-2">
                                        Target Name
                                    </label>
                                    <Input
                                        placeholder="예: HER2, TROP-2, LIV-1..."
                                        value={targetName}
                                        onChange={(e) => setField('targetName', e.target.value)}
                                    />
                                </div>

                                <div>
                                    <label className="block text-sm font-medium text-gray-700 mb-2">
                                        Antibody Selection
                                    </label>
                                    <Select
                                        value={antibodyType}
                                        onValueChange={(value) => setField('antibodyType', value)}
                                    >
                                        <SelectTrigger>
                                            <SelectValue placeholder="Select an antibody..." />
                                        </SelectTrigger>
                                        <SelectContent>
                                            {antibodyOptions.map((option) => (
                                                <SelectItem key={option.value} value={option.value}>
                                                    {option.label}
                                                </SelectItem>
                                            ))}
                                        </SelectContent>
                                    </Select>
                                </div>

                                {antibodyType === 'custom' && (
                                    <div>
                                        <label className="block text-sm font-medium text-gray-700 mb-2">
                                            FASTA Sequence
                                        </label>
                                        <Textarea
                                            placeholder={`>Antibody_Heavy_Chain\nMVLQDRSMG...`}
                                            rows={6}
                                            value={customSequence}
                                            onChange={(e) => validateSequence(e.target.value)}
                                            className={sequenceError ? 'border-red-500 focus:ring-red-500' : ''}
                                        />

                                        {/* 실시간 유효성 검사 에러 표시 */}
                                        {sequenceError && (
                                            <motion.div
                                                initial={{ opacity: 0, y: -10 }}
                                                animate={{ opacity: 1, y: 0 }}
                                                className="mt-2 p-3 bg-red-50 border border-red-200 rounded-lg flex items-start gap-2"
                                            >
                                                <AlertCircle className="w-5 h-5 text-red-500 flex-shrink-0 mt-0.5" />
                                                <div>
                                                    <p className="text-sm font-medium text-red-700">{sequenceError}</p>
                                                    <p className="text-xs text-red-500 mt-1">
                                                        유효한 Amino Acid 코드: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
                                                    </p>
                                                </div>
                                            </motion.div>
                                        )}

                                        {!sequenceError && customSequence && (
                                            <motion.div
                                                initial={{ opacity: 0 }}
                                                animate={{ opacity: 1 }}
                                                className="mt-2 p-2 bg-green-50 border border-green-200 rounded-lg flex items-center gap-2"
                                            >
                                                <Check className="w-4 h-4 text-green-500" />
                                                <span className="text-sm text-green-700">유효한 서열입니다</span>
                                            </motion.div>
                                        )}

                                        <div className="mt-3">
                                            <Button variant="outline" size="sm">
                                                <Upload className="w-4 h-4 mr-2" />
                                                .pdb / .fasta 파일 업로드
                                            </Button>
                                        </div>
                                    </div>
                                )}
                            </CardContent>
                        </Card>
                    )}

                    {currentStep === 2 && (
                        <Card>
                            <CardHeader>
                                <CardTitle>Step 2: Payload & Linker</CardTitle>
                                <CardDescription>약물(Payload)과 링커를 선택하세요</CardDescription>
                            </CardHeader>
                            <CardContent className="space-y-6">
                                <div>
                                    <label className="block text-sm font-medium text-gray-700 mb-2">
                                        Payload Selection
                                    </label>
                                    {/* 2D 화학 구조 미리보기 포함 카드 */}
                                    <div className="grid grid-cols-3 gap-4">
                                        {payloadOptions.map((option) => (
                                            <div
                                                key={option.value}
                                                onClick={() => setField('payloadId', option.value)}
                                                className={`p-4 rounded-xl border-2 cursor-pointer transition-all ${payloadId === option.value
                                                    ? 'border-[#007AFF] bg-blue-50 shadow-md'
                                                    : 'border-gray-200 hover:border-gray-300 hover:shadow-sm'
                                                    }`}
                                            >
                                                {/* 2D 화학 구조 SVG */}
                                                <div
                                                    className="mb-3 p-2 bg-white rounded-lg border"
                                                    dangerouslySetInnerHTML={{ __html: option.structure }}
                                                />
                                                <p className="font-semibold text-gray-900">{option.label}</p>
                                                <p className="text-xs text-gray-500">{option.fullName}</p>
                                                <div className="mt-2 flex items-center justify-between">
                                                    <Badge variant="outline" className="text-xs">
                                                        {option.category.split(' ')[0]}
                                                    </Badge>
                                                    <span className="text-xs text-gray-400">MW: {option.mw}</span>
                                                </div>
                                            </div>
                                        ))}
                                    </div>
                                </div>

                                <div>
                                    <label className="block text-sm font-medium text-gray-700 mb-2">
                                        Linker Selection
                                    </label>
                                    <div className="grid grid-cols-3 gap-3">
                                        {linkerOptions.map((option) => (
                                            <div
                                                key={option.value}
                                                onClick={() => setField('linkerId', option.value)}
                                                className={`p-3 rounded-lg border-2 cursor-pointer transition-all ${linkerId === option.value
                                                    ? 'border-[#007AFF] bg-blue-50'
                                                    : 'border-gray-200 hover:border-gray-300'
                                                    }`}
                                            >
                                                <p className="font-medium text-gray-900 text-sm">{option.label}</p>
                                                <p className="text-xs text-gray-500">{option.desc}</p>
                                            </div>
                                        ))}
                                    </div>
                                </div>

                                <div>
                                    <label className="block text-sm font-medium text-gray-700 mb-2">
                                        DAR (Drug-to-Antibody Ratio): <span className="font-bold text-[#007AFF]">{dar}</span>
                                    </label>
                                    <input
                                        type="range"
                                        min={1}
                                        max={8}
                                        value={dar}
                                        onChange={(e) => setField('dar', parseInt(e.target.value))}
                                        className="w-full h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer accent-[#007AFF]"
                                    />
                                    <div className="flex justify-between text-xs text-gray-500 mt-1">
                                        <span>1</span>
                                        <span className="text-[#007AFF] font-medium">4 (Recommended)</span>
                                        <span>8</span>
                                    </div>
                                </div>
                            </CardContent>
                        </Card>
                    )}

                    {currentStep === 3 && (
                        <Card>
                            <CardHeader>
                                <CardTitle>Step 3: Configuration</CardTitle>
                                <CardDescription>시뮬레이션 설정을 확인하고 실행하세요</CardDescription>
                            </CardHeader>
                            <CardContent className="space-y-6">
                                <div>
                                    <label className="block text-sm font-medium text-gray-700 mb-2">
                                        Simulation Mode
                                    </label>
                                    <div className="grid grid-cols-2 gap-4">
                                        <div
                                            onClick={() => setField('mode', 'fast')}
                                            className={`p-4 rounded-lg border-2 cursor-pointer transition-all ${mode === 'fast'
                                                ? 'border-[#007AFF] bg-blue-50'
                                                : 'border-gray-200 hover:border-gray-300'
                                                }`}
                                        >
                                            <div className="flex items-center justify-between mb-2">
                                                <span className="font-medium text-gray-900">Fast Scan</span>
                                                <Badge variant="outline">1 Credit</Badge>
                                            </div>
                                            <p className="text-sm text-gray-500">3D 구조만 빠르게 확인</p>
                                        </div>
                                        <div
                                            onClick={() => setField('mode', 'deep')}
                                            className={`p-4 rounded-lg border-2 cursor-pointer transition-all ${mode === 'deep'
                                                ? 'border-[#007AFF] bg-blue-50'
                                                : 'border-gray-200 hover:border-gray-300'
                                                }`}
                                        >
                                            <div className="flex items-center justify-between mb-2">
                                                <span className="font-medium text-gray-900">Deep Analysis</span>
                                                <Badge className="bg-[#007AFF]">10 Credits</Badge>
                                            </div>
                                            <p className="text-sm text-gray-500">독성, 특허, 경쟁사 분석 포함</p>
                                        </div>
                                    </div>
                                </div>

                                <div>
                                    <label className="block text-sm font-medium text-gray-700 mb-2">
                                        Job Name
                                    </label>
                                    <Input
                                        placeholder="예: LIV-1_MMAE_Test_01"
                                        value={jobName}
                                        onChange={(e) => setField('jobName', e.target.value)}
                                    />
                                </div>

                                {/* Summary */}
                                <div className="bg-gray-50 rounded-lg p-4">
                                    <h4 className="font-medium text-gray-900 mb-3">Summary</h4>
                                    <div className="grid grid-cols-2 gap-2 text-sm">
                                        <span className="text-gray-500">Target:</span>
                                        <span className="text-gray-900">{targetName || '-'}</span>
                                        <span className="text-gray-500">Antibody:</span>
                                        <span className="text-gray-900">{antibodyType || '-'}</span>
                                        <span className="text-gray-500">Payload:</span>
                                        <span className="text-gray-900">{payloadId || '-'}</span>
                                        <span className="text-gray-500">Linker:</span>
                                        <span className="text-gray-900">{linkerId || '-'}</span>
                                        <span className="text-gray-500">DAR:</span>
                                        <span className="text-gray-900">{dar}</span>
                                        <span className="text-gray-500">Mode:</span>
                                        <span className="text-gray-900 font-medium">
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
                >
                    <ChevronLeft className="w-4 h-4 mr-2" />
                    Back
                </Button>

                {currentStep < 3 ? (
                    <Button
                        onClick={handleNext}
                        className="bg-[#007AFF] hover:bg-[#0056b3]"
                        disabled={currentStep === 1 && antibodyType === 'custom' && !!sequenceError}
                    >
                        Next
                        <ChevronRight className="w-4 h-4 ml-2" />
                    </Button>
                ) : (
                    <Button onClick={handleSubmit} className="bg-[#007AFF] hover:bg-[#0056b3]">
                        <Rocket className="w-4 h-4 mr-2" />
                        Run Simulation
                    </Button>
                )}
            </div>
        </div>
    )
}
