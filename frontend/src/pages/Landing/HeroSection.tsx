import { motion } from 'framer-motion'
import { Button } from '@/components/ui/button'
import { Link } from 'react-router-dom'
import { ArrowRight, Sparkles } from 'lucide-react'

export function HeroSection() {
    return (
        <section className="relative min-h-screen flex items-center bg-[#0F172A] overflow-hidden pt-20">
            <div className="w-full max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
                <div className="grid grid-cols-1 lg:grid-cols-2 gap-12 items-center">
                    {/* Left Column: Text Content */}
                    <motion.div
                        initial={{ opacity: 0, x: -20 }}
                        animate={{ opacity: 1, x: 0 }}
                        transition={{ duration: 0.8 }}
                        className="text-left"
                    >
                        <div className="inline-flex items-center gap-2 px-4 py-2 bg-blue-500/10 rounded-full mb-8 border border-blue-500/20">
                            <Sparkles className="w-4 h-4 text-[#007AFF]" />
                            <span className="text-sm text-blue-400 font-medium tracking-wide">AstraForge: The Brain of ADC Discovery</span>
                        </div>

                        <h1 className="text-4xl md:text-6xl font-bold text-white leading-tight mb-6 tracking-tight">
                            Accelerate ADC Discovery <br />
                            <span className="text-transparent bg-clip-text bg-gradient-to-r from-[#007AFF] to-[#06B6D4]">
                                from Months to Minutes
                            </span>
                        </h1>

                        <p className="text-lg text-gray-400 mb-8 max-w-xl leading-relaxed">
                            실험실의 시행착오를 AI 시뮬레이션으로 대체하십시오. <br />
                            AstraForge는 FDA 승인 데이터를 학습한 6개의 AI 에이전트로 귀사의 링커-페이로드 최적화를 가속화합니다.
                        </p>

                        <div className="flex flex-col sm:flex-row gap-4 mb-12">
                            <Button size="lg" className="bg-[#007AFF] hover:bg-[#0056b3] text-white px-8 py-6 text-lg rounded-xl" asChild>
                                <Link to="/signup">
                                    Start Research
                                    <ArrowRight className="w-5 h-5 ml-2" />
                                </Link>
                            </Button>
                            <Button size="lg" variant="outline" className="text-white border-white/10 hover:bg-white/5 px-8 py-6 text-lg rounded-xl" asChild>
                                <Link to="/demo">Request Demo</Link>
                            </Button>
                        </div>

                        {/* Social Proof */}
                        <div className="pt-8 border-t border-white/10">
                            <p className="text-xs text-gray-500 mb-4 font-semibold tracking-widest uppercase">TRUSTED DATA SOURCES</p>
                            <div className="flex gap-6 opacity-50 grayscale">
                                <span className="text-white font-bold text-lg">FDA</span>
                                <span className="text-white font-bold text-lg">PubMed</span>
                                <span className="text-white font-bold text-lg">NIH</span>
                                <span className="text-white font-bold text-lg">USPTO</span>
                            </div>
                        </div>
                    </motion.div>

                    {/* Right Column: Video */}
                    <motion.div
                        initial={{ opacity: 0, x: 20 }}
                        animate={{ opacity: 1, x: 0 }}
                        transition={{ duration: 0.8, delay: 0.2 }}
                        className="relative"
                    >
                        <div className="relative rounded-2xl overflow-hidden shadow-2xl shadow-blue-900/20 border border-white/10 bg-slate-900">
                            <video
                                autoPlay
                                loop
                                muted
                                playsInline
                                className="w-full h-auto"
                            >
                                <source src="/1.mp4" type="video/mp4" />
                            </video>
                            {/* Optional Overlay for better text contrast if needed, but usually not for side video */}
                            <div className="absolute inset-0 bg-gradient-to-tr from-black/20 to-transparent pointer-events-none" />
                        </div>
                        {/* Decorative Elements */}
                        <div className="absolute -z-10 top-1/2 left-1/2 -translate-x-1/2 -translate-y-1/2 w-[120%] h-[120%] bg-blue-500/20 blur-[100px] rounded-full opacity-30" />
                    </motion.div>
                </div>
            </div>
        </section>
    )
}
