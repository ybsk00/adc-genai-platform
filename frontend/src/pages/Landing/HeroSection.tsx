import { motion } from 'framer-motion'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { ArrowRight, Sparkles } from 'lucide-react'
import { useState } from 'react'

export function HeroSection() {
    const [email, setEmail] = useState('')

    const handleSubmit = (e: React.FormEvent) => {
        e.preventDefault()
        // TODO: API 호출 - 골든셋 PDF 다운로드 링크 발송
        console.log('Lead magnet email:', email)
    }

    return (
        <section className="relative min-h-screen flex items-center pt-16 overflow-hidden bg-gradient-to-br from-[#0A0F1C] via-[#1A1F2E] to-[#0A0F1C]">
            {/* Background Effects */}
            <div className="absolute inset-0 overflow-hidden">
                <div className="absolute -top-40 -right-40 w-80 h-80 bg-[#007AFF]/20 rounded-full blur-3xl" />
                <div className="absolute -bottom-40 -left-40 w-80 h-80 bg-purple-500/20 rounded-full blur-3xl" />
            </div>

            <div className="relative max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-20">
                <div className="grid lg:grid-cols-2 gap-12 items-center">
                    {/* Left: Text Content */}
                    <motion.div
                        initial={{ opacity: 0, x: -50 }}
                        animate={{ opacity: 1, x: 0 }}
                        transition={{ duration: 0.8 }}
                    >
                        <div className="inline-flex items-center gap-2 px-4 py-2 bg-white/10 rounded-full mb-6">
                            <Sparkles className="w-4 h-4 text-[#007AFF]" />
                            <span className="text-sm text-gray-300">AI-Powered Drug Discovery</span>
                        </div>

                        <h1 className="text-4xl md:text-5xl lg:text-6xl font-bold text-white leading-tight mb-6">
                            ADC 개발의 성공,
                            <br />
                            <span className="text-transparent bg-clip-text bg-gradient-to-r from-[#007AFF] to-purple-500">
                                '실험'하지 마시고
                            </span>
                            <br />
                            '설계'하세요.
                        </h1>

                        <p className="text-lg text-gray-400 mb-8 max-w-lg">
                            AI 기반 시뮬레이션으로 Linker-Payload 최적화 및 독성 예측을
                            몇 달에서 몇 분으로 단축하세요.
                        </p>

                        {/* Email Capture Form */}
                        <form onSubmit={handleSubmit} className="flex flex-col sm:flex-row gap-3 max-w-md">
                            <Input
                                type="email"
                                placeholder="업무용 이메일을 입력하세요"
                                value={email}
                                onChange={(e) => setEmail(e.target.value)}
                                className="bg-white/10 border-white/20 text-white placeholder:text-gray-500"
                            />
                            <Button
                                type="submit"
                                className="bg-[#007AFF] hover:bg-[#0056b3] whitespace-nowrap"
                            >
                                Golden Set 받기
                                <ArrowRight className="w-4 h-4 ml-2" />
                            </Button>
                        </form>

                        <p className="text-sm text-gray-500 mt-3">
                            무료로 FDA 승인 ADC 15종 분석 리포트를 받아보세요.
                        </p>
                    </motion.div>

                    {/* Right: 3D Visualization Placeholder */}
                    <motion.div
                        initial={{ opacity: 0, x: 50 }}
                        animate={{ opacity: 1, x: 0 }}
                        transition={{ duration: 0.8, delay: 0.2 }}
                        className="relative"
                    >
                        <div className="aspect-square max-w-lg mx-auto relative">
                            {/* Placeholder for MolStar 3D Viewer */}
                            <div className="absolute inset-0 bg-gradient-to-br from-[#007AFF]/20 to-purple-500/20 rounded-3xl backdrop-blur-sm border border-white/10">
                                <div className="absolute inset-0 flex items-center justify-center">
                                    <div className="text-center">
                                        <div className="w-32 h-32 mx-auto mb-4 rounded-full bg-gradient-to-br from-[#007AFF] to-purple-500 animate-pulse flex items-center justify-center">
                                            <span className="text-4xl">🧬</span>
                                        </div>
                                        <p className="text-gray-400 text-sm">Interactive 3D ADC Model</p>
                                        <p className="text-gray-500 text-xs mt-1">(MolStar Integration)</p>
                                    </div>
                                </div>
                            </div>

                            {/* Floating Stats */}
                            <motion.div
                                animate={{ y: [0, -10, 0] }}
                                transition={{ duration: 3, repeat: Infinity }}
                                className="absolute -top-4 -right-4 px-4 py-2 bg-white/10 backdrop-blur-md rounded-lg border border-white/20"
                            >
                                <p className="text-xs text-gray-400">분석 정확도</p>
                                <p className="text-lg font-bold text-white">95.7%</p>
                            </motion.div>

                            <motion.div
                                animate={{ y: [0, 10, 0] }}
                                transition={{ duration: 3, repeat: Infinity, delay: 0.5 }}
                                className="absolute -bottom-4 -left-4 px-4 py-2 bg-white/10 backdrop-blur-md rounded-lg border border-white/20"
                            >
                                <p className="text-xs text-gray-400">분석 시간</p>
                                <p className="text-lg font-bold text-white">&lt; 5분</p>
                            </motion.div>
                        </div>
                    </motion.div>
                </div>
            </div>
        </section>
    )
}
