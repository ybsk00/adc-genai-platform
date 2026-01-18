import { motion } from 'framer-motion'
import { Button } from '@/components/ui/button'
import { ArrowRight, Sparkles } from 'lucide-react'

export function HeroSection() {
    return (
        <section className="relative min-h-screen flex items-center pt-16 overflow-hidden bg-[#0F172A]">
            {/* Background Video with Overlay */}
            <div className="absolute inset-0 overflow-hidden">
                <video
                    autoPlay
                    loop
                    muted
                    playsInline
                    className="absolute w-full h-full object-cover"
                >
                    <source src="/1.mp4" type="video/mp4" />
                </video>
                {/* Black Overlay (Opacity 50%) */}
                <div className="absolute inset-0 bg-black/50" />
            </div>

            <div className="relative max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-20">
                <div className="grid lg:grid-cols-2 gap-12 items-center">
                    {/* Left: Text Content */}
                    <motion.div
                        initial={{ opacity: 0, x: -50 }}
                        animate={{ opacity: 1, x: 0 }}
                        transition={{ duration: 0.8 }}
                    >
                        <div className="inline-flex items-center gap-2 px-4 py-2 bg-white/10 rounded-full mb-6 backdrop-blur-sm border border-white/10">
                            <Sparkles className="w-4 h-4 text-[#06B6D4]" />
                            <span className="text-sm text-gray-200">AstraForge: The Brain of ADC Discovery</span>
                        </div>

                        <h1 className="text-4xl md:text-5xl lg:text-6xl font-bold text-white leading-tight mb-6">
                            Accelerate ADC Discovery
                            <br />
                            <span className="text-transparent bg-clip-text bg-gradient-to-r from-[#06B6D4] to-[#10B981]">
                                from Months to Minutes.
                            </span>
                        </h1>

                        <p className="text-lg text-gray-200 mb-8 max-w-lg font-light">
                            실험실의 시행착오를 AI 시뮬레이션으로 대체하십시오.
                            <br />
                            AstraForge는 FDA 승인 데이터를 학습한 6개의 AI 에이전트로 귀사의 링커-페이로드 최적화를 가속화합니다.
                        </p>

                        <div className="flex flex-col sm:flex-row gap-4">
                            <Button
                                size="lg"
                                className="bg-[#06B6D4] hover:bg-[#0891B2] text-white font-semibold px-8"
                            >
                                Start Research
                                <ArrowRight className="w-4 h-4 ml-2" />
                            </Button>
                            <Button
                                size="lg"
                                variant="outline"
                                className="border-white/20 text-white hover:bg-white/10"
                            >
                                Request Demo
                            </Button>
                        </div>

                        {/* Social Proof */}
                        <div className="mt-12 pt-8 border-t border-white/10">
                            <p className="text-sm text-gray-400 mb-4 font-semibold tracking-wider">TRUSTED DATA</p>
                            <div className="flex flex-wrap gap-6 items-center text-gray-500 text-sm font-mono">
                                <span>Trained on 2,000+ Clinical Trials</span>
                                <span className="w-1 h-1 bg-gray-600 rounded-full" />
                                <span>15+ FDA Approved ADCs</span>
                            </div>
                            <div className="flex gap-6 mt-4 opacity-50 grayscale">
                                {/* Text placeholders for logos as requested */}
                                <span className="text-white font-bold text-lg">FDA</span>
                                <span className="text-white font-bold text-lg">PubMed</span>
                                <span className="text-white font-bold text-lg">NIH</span>
                                <span className="text-white font-bold text-lg">USPTO</span>
                            </div>
                        </div>
                    </motion.div>

                    {/* Right: Empty for now as video is background, or we can keep the 3D placeholder if desired, 
                        but the request implies the video is the main visual. 
                        Let's keep the layout balanced but maybe remove the placeholder if it conflicts with the video.
                        However, the user asked for "Video Background", so the right side might be better left empty or used for a glass card.
                        For now, I'll remove the placeholder to let the video shine, or maybe add a small glass card stats.
                    */}
                    <motion.div
                        initial={{ opacity: 0, x: 50 }}
                        animate={{ opacity: 1, x: 0 }}
                        transition={{ duration: 0.8, delay: 0.2 }}
                        className="hidden lg:block relative"
                    >
                        {/* Optional: Glass Card Stats floating */}
                        <div className="absolute right-0 top-1/2 -translate-y-1/2 w-72 p-6 bg-black/30 backdrop-blur-md rounded-2xl border border-white/10">
                            <div className="space-y-4">
                                <div>
                                    <p className="text-xs text-gray-400 uppercase tracking-wider">Analysis Speed</p>
                                    <p className="text-2xl font-bold text-[#06B6D4]">100x Faster</p>
                                </div>
                                <div>
                                    <p className="text-xs text-gray-400 uppercase tracking-wider">Accuracy</p>
                                    <p className="text-2xl font-bold text-[#10B981]">98.5%</p>
                                </div>
                                <div>
                                    <p className="text-xs text-gray-400 uppercase tracking-wider">Cost Reduction</p>
                                    <p className="text-2xl font-bold text-purple-400">~80%</p>
                                </div>
                            </div>
                        </div>
                    </motion.div>
                </div>
            </div>
        </section>
    )
}
