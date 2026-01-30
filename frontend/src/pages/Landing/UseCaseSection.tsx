import { motion } from 'framer-motion'
import ForBioResearchers from '@/assets/icons/For_Bio-Researchers.png'
import ForVCInvestors from '@/assets/icons/For_VC_Investors.png'
import ForCDMOs from '@/assets/icons/For_CDMOs.png'

export function UseCaseSection() {
    return (
        <section className="py-24 bg-[#0F172A] relative">
            <div className="text-center mb-16">
                <h2 className="text-3xl md:text-5xl font-bold text-white mb-6">
                    Solutions by Role
                </h2>
                <p className="text-gray-400 text-lg max-w-2xl mx-auto">
                    Tailored solutions for experts in every field.
                </p>
            </div>

            <div className="grid grid-cols-1 md:grid-cols-3 gap-8 max-w-7xl mx-auto px-4">
                {/* Bio-Researcher */}
                <div className="group relative p-8 rounded-3xl bg-[#0F172A]/70 backdrop-blur-[20px] border border-white/10 hover:border-blue-500/50 transition-all duration-300 hover:-translate-y-2">
                    <div className="absolute inset-0 bg-gradient-to-b from-blue-500/5 to-transparent opacity-0 group-hover:opacity-100 transition-opacity duration-300 rounded-3xl" />

                    <div className="relative z-10">
                        <div className="w-14 h-14 rounded-2xl bg-blue-500/10 flex items-center justify-center mb-6 group-hover:bg-blue-500/20 transition-colors overflow-hidden">
                            <img src={ForBioResearchers} alt="For Bio-Researchers" className="w-full h-full object-cover" />
                        </div>

                        <h3 className="text-xl font-bold text-white mb-2">For Bio-Researchers</h3>
                        <p className="text-blue-400 text-sm font-medium mb-4">Researchers</p>

                        <p className="text-gray-400 leading-relaxed">
                            "Validate toxicity before synthesis to drastically reduce wet-lab experiment cycles and optimize candidate screening."
                        </p>

                        <ul className="mt-6 space-y-3">
                            <li className="flex items-center text-sm text-gray-500 group-hover:text-gray-300 transition-colors">
                                <span className="w-1.5 h-1.5 rounded-full bg-blue-500 mr-2" />
                                Optimize Experiment Design
                            </li>
                            <li className="flex items-center text-sm text-gray-500 group-hover:text-gray-300 transition-colors">
                                <span className="w-1.5 h-1.5 rounded-full bg-blue-500 mr-2" />
                                Accelerate Candidate Screening
                            </li>
                        </ul>
                    </div>
                </div>

                {/* VC Investor */}
                <div className="group relative p-8 rounded-3xl bg-[#0F172A]/70 backdrop-blur-[20px] border border-white/10 hover:border-purple-500/50 transition-all duration-300 hover:-translate-y-2">
                    <div className="absolute inset-0 bg-gradient-to-b from-purple-500/5 to-transparent opacity-0 group-hover:opacity-100 transition-opacity duration-300 rounded-3xl" />

                    <div className="relative z-10">
                        <div className="w-14 h-14 rounded-2xl bg-purple-500/10 flex items-center justify-center mb-6 group-hover:bg-purple-500/20 transition-colors overflow-hidden">
                            <img src={ForVCInvestors} alt="For VC Investors" className="w-full h-full object-cover" />
                        </div>

                        <h3 className="text-xl font-bold text-white mb-2">For VC Investors</h3>
                        <p className="text-purple-400 text-sm font-medium mb-4">Venture Capital</p>

                        <p className="text-gray-400 leading-relaxed">
                            "Verify technical feasibility and patent risks in minutes to make data-driven, high-confidence investment decisions."
                        </p>

                        <ul className="mt-6 space-y-3">
                            <li className="flex items-center text-sm text-gray-500 group-hover:text-gray-300 transition-colors">
                                <span className="w-1.5 h-1.5 rounded-full bg-purple-500 mr-2" />
                                Technical Valuation
                            </li>
                            <li className="flex items-center text-sm text-gray-500 group-hover:text-gray-300 transition-colors">
                                <span className="w-1.5 h-1.5 rounded-full bg-purple-500 mr-2" />
                                FTO Risk Analysis
                            </li>
                        </ul>
                    </div>
                </div>

                {/* CDMO */}
                <div className="group relative p-8 rounded-3xl bg-[#0F172A]/70 backdrop-blur-[20px] border border-white/10 hover:border-green-500/50 transition-all duration-300 hover:-translate-y-2">
                    <div className="absolute inset-0 bg-gradient-to-b from-green-500/5 to-transparent opacity-0 group-hover:opacity-100 transition-opacity duration-300 rounded-3xl" />

                    <div className="relative z-10">
                        <div className="w-14 h-14 rounded-2xl bg-green-500/10 flex items-center justify-center mb-6 group-hover:bg-green-500/20 transition-colors overflow-hidden">
                            <img src={ForCDMOs} alt="For CDMOs" className="w-full h-full object-cover" />
                        </div>

                        <h3 className="text-xl font-bold text-white mb-2">For CDMOs</h3>
                        <p className="text-green-400 text-sm font-medium mb-4">Contract Development & Manufacturing</p>

                        <p className="text-gray-400 leading-relaxed">
                            "Propose optimal linker/payload combinations with verified reasoning to increase client trust and development efficiency."
                        </p>

                        <ul className="mt-6 space-y-3">
                            <li className="flex items-center text-sm text-gray-500 group-hover:text-gray-300 transition-colors">
                                <span className="w-1.5 h-1.5 rounded-full bg-green-500 mr-2" />
                                Customized Proposals
                            </li>
                            <li className="flex items-center text-sm text-gray-500 group-hover:text-gray-300 transition-colors">
                                <span className="w-1.5 h-1.5 rounded-full bg-green-500 mr-2" />
                                Increase Development Efficiency
                            </li>
                        </ul>
                    </div>
                </div>
            </div>
        </section>
    )
}
