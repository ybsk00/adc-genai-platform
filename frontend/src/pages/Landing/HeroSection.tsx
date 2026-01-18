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

            <div className="relative max-w-5xl mx-auto px-4 sm:px-6 lg:px-8 py-32 md:py-40 text-center">
                <motion.div
                    initial={{ opacity: 0, y: 30 }}
                    animate={{ opacity: 1, y: 0 }}
                    transition={{ duration: 0.8 }}
                    className="flex flex-col items-center"
                >
                    <div className="inline-flex items-center gap-2 px-4 py-2 bg-white/10 rounded-full mb-8 backdrop-blur-sm border border-white/10">
                        <Sparkles className="w-4 h-4 text-[#06B6D4]" />
                        <span className="text-sm text-gray-200 tracking-wide">AstraForge: The Brain of ADC Discovery</span>
                    </div>

                    <h1 className="text-4xl md:text-6xl lg:text-7xl font-bold text-white leading-tight mb-8 tracking-tight">
                        Accelerate ADC Discovery
                        <br className="hidden md:block" />
                        <span className="text-transparent bg-clip-text bg-gradient-to-r from-[#06B6D4] to-[#10B981] ml-2">
                            from Months to Minutes.
                        </span>
                    </h1>

                    <p className="text-lg md:text-xl text-gray-200 mb-10 max-w-2xl font-light leading-relaxed tracking-wide">
                        실험실의 시행착오를 AI 시뮬레이션으로 대체하십시오.
                        <br />
                        AstraForge는 FDA 승인 데이터를 학습한 6개의 AI 에이전트로 귀사의 링커-페이로드 최적화를 가속화합니다.
                    </p>

                    <div className="flex flex-col sm:flex-row gap-6 mb-16">
                        <Button
                            size="lg"
                            className="bg-[#06B6D4] hover:bg-[#0891B2] text-white font-semibold px-8 py-6 text-lg"
                        >
                            Start Research
                            <ArrowRight className="w-5 h-5 ml-2" />
                        </Button>
                        <Button
                            size="lg"
                            variant="outline"
                            className="border-white/20 text-white hover:bg-white/10 px-8 py-6 text-lg"
                        >
                            Request Demo
                        </Button>
                    </div>

                    {/* Social Proof */}
                    <div className="pt-10 border-t border-white/10 w-full max-w-3xl">
                        <p className="text-sm text-gray-400 mb-6 font-semibold tracking-widest uppercase">TRUSTED DATA</p>
                        <div className="flex flex-wrap justify-center gap-8 items-center text-gray-500 text-sm font-mono mb-6">
                            <span>Trained on 2,000+ Clinical Trials</span>
                            <span className="w-1.5 h-1.5 bg-gray-600 rounded-full" />
                            <span>15+ FDA Approved ADCs</span>
                        </div>
                        <div className="flex justify-center gap-10 opacity-50 grayscale">
                            <span className="text-white font-bold text-xl tracking-wider">FDA</span>
                            <span className="text-white font-bold text-xl tracking-wider">PubMed</span>
                            <span className="text-white font-bold text-xl tracking-wider">NIH</span>
                            <span className="text-white font-bold text-xl tracking-wider">USPTO</span>
                        </div>
                    </div>
                </motion.div>
            </div>
        </section>
    )
}
