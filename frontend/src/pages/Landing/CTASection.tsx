import { Button } from '@/components/ui/button'
import { Link } from 'react-router-dom'
import { ArrowRight, Download } from 'lucide-react'

export function CTASection() {
    return (
        <section className="py-32 relative overflow-hidden">
            {/* Background Gradient */}
            <div className="absolute inset-0 bg-gradient-to-b from-[#0F172A] to-blue-900/20" />
            <div className="absolute inset-0 bg-[url('/grid.svg')] opacity-20" />

            <div className="relative z-10 max-w-4xl mx-auto px-4 text-center">
                <h2 className="text-4xl md:text-6xl font-bold text-white mb-6 tracking-tight">
                    Ready to Transform Your <br />
                    <span className="text-transparent bg-clip-text bg-gradient-to-r from-blue-400 to-cyan-400">
                        ADC Research?
                    </span>
                </h2>

                <p className="text-xl text-gray-400 mb-12 leading-relaxed">
                    Request an enterprise demo now and receive <br className="hidden md:block" />
                    a customized analysis report for your target.
                </p>

                <div className="flex flex-col sm:flex-row gap-4 justify-center">
                    <Button size="lg" className="bg-[#007AFF] hover:bg-[#0056b3] text-white px-8 h-14 rounded-xl text-lg font-semibold shadow-lg shadow-blue-500/20 transition-all hover:scale-105" asChild>
                        <Link to="/signup">
                            Request Enterprise Demo <ArrowRight className="ml-2 w-5 h-5" />
                        </Link>
                    </Button>
                    <Button size="lg" variant="outline" className="border-white/10 text-white hover:bg-white/5 h-14 rounded-xl text-lg font-medium px-8" asChild>
                        <Link to="/report-sample">
                            Download Report Sample <Download className="ml-2 w-5 h-5" />
                        </Link>
                    </Button>
                </div>
            </div>
        </section>
    )
}
