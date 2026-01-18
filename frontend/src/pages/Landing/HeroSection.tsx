import { motion } from 'framer-motion'
import { Button } from '@/components/ui/button'
import { Link } from 'react-router-dom'
import { Container } from '@/components/layout/Container'

export function HeroSection() {
    return (
        <section className="relative min-h-screen flex items-center justify-center pt-16 overflow-hidden bg-[#0F172A]" style={{ display: 'flex', justifyContent: 'center', alignItems: 'center' }}>
            {/* Background Video with Overlay */}
            <div className="absolute inset-0 overflow-hidden">
                <video
                    autoPlay
                    loop
                    muted
                    playsInline
                    className="w-full h-full object-cover opacity-30"
                >
                    <source src="/1.mp4" type="video/mp4" />
                </video>
                <div className="absolute inset-0 bg-black/50" />
            </div>

            <Container className="relative py-32 md:py-40 text-center">
                <motion.div
                    initial={{ opacity: 0, y: 30 }}
                    animate={{ opacity: 1, y: 0 }}
                    transition={{ duration: 0.8 }}
                    className="max-w-4xl mx-auto"
                >
                    <h1 className="text-5xl md:text-7xl font-bold text-white mb-8 tracking-tight leading-tight">
                        Revolutionizing ADC Discovery <br />
                        <span className="text-transparent bg-clip-text bg-gradient-to-r from-blue-400 to-cyan-300">
                            with Generative AI
                        </span>
                    </h1>
                    <p className="text-xl md:text-2xl text-gray-300 mb-12 max-w-2xl mx-auto leading-relaxed">
                        From Target Identification to Toxicity Prediction. <br />
                        Accelerate your research from months to minutes.
                    </p>
                    <div className="flex flex-col sm:flex-row gap-4 justify-center items-center">
                        <Button size="lg" className="bg-[#007AFF] hover:bg-[#0056b3] text-lg px-8 py-6 h-auto rounded-full shadow-lg shadow-blue-500/25 transition-all hover:scale-105" asChild>
                            <Link to="/signup">Start Free Trial</Link>
                        </Button>
                        <Button size="lg" variant="outline" className="text-white border-white/20 hover:bg-white/10 text-lg px-8 py-6 h-auto rounded-full backdrop-blur-sm transition-all hover:scale-105" asChild>
                            <Link to="/demo">Watch Demo</Link>
                        </Button>
                    </div>
                </motion.div>
            </Container>
        </section>
    )
}
