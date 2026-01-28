import { Button } from '@/components/ui/button'
import { Link } from 'react-router-dom'
import { ArrowRight, Download } from 'lucide-react'

export function HeroSection() {
    return (
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-12 items-center">
            {/* Left: Text Content */}
            <div className="space-y-8 z-10">

                {/* Stats Ticker */}
                <div className="inline-flex items-center gap-3 px-4 py-2 rounded-full bg-blue-500/10 border border-blue-500/20 text-blue-400 text-xs font-medium animate-fade-in-up">
                    <span className="flex items-center gap-1"><span className="w-2 h-2 rounded-full bg-green-500 animate-pulse" /> Live: 15,420 Papers Indexed</span>
                    <span className="w-px h-3 bg-blue-500/30" />
                    <span className="flex items-center gap-1">🧬 152 FDA Approved ADCs</span>
                    <span className="w-px h-3 bg-blue-500/30" />
                    <span className="flex items-center gap-1">⚡ 99.8% Prediction Accuracy</span>
                </div>

                <h1 className="text-5xl lg:text-7xl font-bold tracking-tight leading-tight text-white">
                    Accelerate ADC Discovery: <br />
                    From "Black Box" to <br />
                    <span className="text-transparent bg-clip-text bg-gradient-to-r from-blue-400 to-cyan-400">
                        "Verified Reasoning"
                    </span>
                </h1>

                <p className="text-lg text-gray-400 max-w-xl leading-relaxed">
                    Replace lab trial and error with Hallucination-Free AI Simulation. <br />
                    AstraForge ensures every design is physically valid through real-time Python/RDKit code execution.
                </p>

                <div className="flex flex-wrap gap-4">
                    <Button size="lg" className="bg-[#007AFF] hover:bg-[#0056b3] text-white px-8 h-12 rounded-xl text-lg font-semibold shadow-lg shadow-blue-500/20 transition-all hover:scale-105" asChild>
                        <Link to="/signup">
                            Start Research <ArrowRight className="ml-2 w-5 h-5" />
                        </Link>
                    </Button>
                    <Button size="lg" variant="outline" className="border-white/10 text-white hover:bg-white/5 h-12 rounded-xl text-lg font-medium px-8" asChild>
                        <Link to="/signup">Request Demo</Link>
                    </Button>
                    <Button
                        size="lg"
                        variant="outline"
                        className="border-blue-400/30 text-blue-400 hover:bg-blue-400/10 h-12 rounded-xl text-lg font-medium px-8"
                        onClick={() => window.open('#', '_blank')}
                    >
                        <Download className="mr-2 w-5 h-5" />
                        Sample Report
                    </Button>
                </div>

                <div className="pt-8 border-t border-white/10">
                    <p className="text-xs font-bold text-gray-500 tracking-wider uppercase mb-4">Trusted Data Sources</p>
                    <div className="flex gap-6 opacity-50 grayscale hover:grayscale-0 transition-all duration-500">
                        <span className="text-lg font-bold text-white">FDA</span>
                        <span className="text-lg font-bold text-white">PubMed</span>
                        <span className="text-lg font-bold text-white">NIH</span>
                        <span className="text-lg font-bold text-white">USPTO</span>
                    </div>
                </div>
            </div>

            {/* Right: Video/Visual */}
            <div className="relative z-0 lg:scale-[1.15] origin-center transition-transform duration-700">
                <div className="relative rounded-2xl overflow-hidden shadow-2xl border border-white/10 bg-[#0F172A]/50 backdrop-blur-sm group">
                    <div className="absolute inset-0 bg-gradient-to-tr from-blue-500/10 to-purple-500/10 opacity-0 group-hover:opacity-100 transition-opacity duration-500" />

                    {/* Video Background */}
                    <video
                        autoPlay
                        loop
                        muted
                        playsInline
                        className="w-full h-full object-cover opacity-80 mix-blend-screen"
                    >
                        <source src="/1.mp4" type="video/mp4" />
                    </video>

                    {/* Overlay Gradient */}
                    <div className="absolute inset-0 bg-gradient-to-t from-[#0F172A] via-transparent to-transparent opacity-60" />
                </div>

                {/* Decorative Elements */}
                <div className="absolute -top-10 -right-10 w-32 h-32 bg-blue-500/20 rounded-full blur-3xl animate-pulse" />
                <div className="absolute -bottom-10 -left-10 w-32 h-32 bg-purple-500/20 rounded-full blur-3xl animate-pulse delay-700" />
            </div>
        </div>
    )
}
