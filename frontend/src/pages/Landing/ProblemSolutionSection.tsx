import { Zap, DollarSign, FlaskConical, TrendingDown, CheckCircle } from 'lucide-react'

export function ProblemSolutionSection() {
    return (
        <section className="py-24 relative overflow-hidden">
            <div className="text-center mb-16">
                <h2 className="text-3xl md:text-5xl font-bold text-white mb-6">
                    Why AI for ADC?
                </h2>
                <p className="text-gray-400 text-lg max-w-2xl mx-auto">
                    Experience overwhelming efficiency based on data, beyond the limits of traditional drug discovery.
                </p>
            </div>

            <div className="grid grid-cols-1 lg:grid-cols-2 gap-8 max-w-6xl mx-auto px-4">
                {/* Traditional Wet Lab */}
                <div className="bg-slate-900/50 border border-red-500/20 rounded-3xl p-8 relative overflow-hidden group">
                    <div className="absolute top-0 left-0 w-full h-1 bg-gradient-to-r from-red-500/50 to-transparent" />
                    <h3 className="text-2xl font-bold text-gray-300 mb-8 flex items-center gap-3">
                        <FlaskConical className="w-6 h-6 text-red-400" />
                        Traditional Wet Lab
                    </h3>

                    <div className="space-y-6">
                        <div className="flex items-start gap-4 p-4 rounded-xl bg-red-500/5 border border-red-500/10">
                            <Zap className="w-6 h-6 text-red-400 shrink-0" />
                            <div>
                                <div className="text-red-400 font-bold text-lg">12~18 Months</div>
                                <div className="text-gray-400 text-sm">Long time to derive initial candidates</div>
                            </div>
                        </div>
                        <div className="flex items-start gap-4 p-4 rounded-xl bg-red-500/5 border border-red-500/10">
                            <DollarSign className="w-6 h-6 text-red-400 shrink-0" />
                            <div>
                                <div className="text-red-400 font-bold text-lg">Billions in Cost</div>
                                <div className="text-gray-400 text-sm">High costs due to repetitive experiments and labor</div>
                            </div>
                        </div>
                        <div className="flex items-start gap-4 p-4 rounded-xl bg-red-500/5 border border-red-500/10">
                            <FlaskConical className="w-6 h-6 text-red-400 shrink-0" />
                            <div>
                                <div className="text-red-400 font-bold text-lg">Trial & Error</div>
                                <div className="text-gray-400 text-sm">Repeated failures and re-experiments</div>
                            </div>
                        </div>
                    </div>
                </div>

                {/* AstraForge AI */}
                <div className="bg-blue-900/20 border border-blue-500/30 rounded-3xl p-8 relative overflow-hidden group">
                    <div className="absolute top-0 left-0 w-full h-1 bg-gradient-to-r from-blue-500 to-cyan-500" />
                    <div className="absolute inset-0 bg-blue-500/5 group-hover:bg-blue-500/10 transition-colors duration-500" />

                    <h3 className="text-2xl font-bold text-white mb-8 flex items-center gap-3 relative z-10">
                        <Zap className="w-6 h-6 text-blue-400" />
                        AstraForge AI
                    </h3>

                    <div className="space-y-6 relative z-10">
                        <div className="flex items-start gap-4 p-4 rounded-xl bg-blue-500/10 border border-blue-500/20">
                            <Zap className="w-6 h-6 text-blue-400 shrink-0" />
                            <div>
                                <div className="text-blue-400 font-bold text-lg">2 Weeks</div>
                                <div className="text-gray-300 text-sm">Only 2 weeks to derive optimal candidates</div>
                            </div>
                        </div>
                        <div className="flex items-start gap-4 p-4 rounded-xl bg-blue-500/10 border border-blue-500/20">
                            <TrendingDown className="w-6 h-6 text-blue-400 shrink-0" />
                            <div>
                                <div className="text-blue-400 font-bold text-lg">90% Cost Reduction</div>
                                <div className="text-gray-300 text-sm">Drastic cost reduction by reducing unnecessary experiments</div>
                            </div>
                        </div>
                        <div className="flex items-start gap-4 p-4 rounded-xl bg-blue-500/10 border border-blue-500/20">
                            <CheckCircle className="w-6 h-6 text-blue-400 shrink-0" />
                            <div>
                                <div className="text-blue-400 font-bold text-lg">In-Silico Pre-validation</div>
                                <div className="text-gray-300 text-sm">Maximize success rate with high-accuracy simulation</div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </section>
    )
}
