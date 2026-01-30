import { motion } from 'framer-motion'
import TheAuditor from '@/assets/icons/The_Auditor.png'
import PredictiveToxicology from '@/assets/icons/Predictive_Toxicology.png'
import IPRiskReasoning from '@/assets/icons/IP_Risk_Reasoning.png'
import CompetitiveIntelligence from '@/assets/icons/Competitive_Intelligence.png'
import TheLibrarian from '@/assets/icons/The_Librarian.png'
import InSilico3DModeling from '@/assets/icons/In_Silico_3D_Modeling.png'

const features = [
    {
        icon: TheAuditor,
        title: "The Auditor",
        description: "Eliminates AI hallucinations by converting designs into executable Python code. Validates MW, stability, and synthesis scores in real-time.",
        color: "from-blue-500 to-cyan-500"
    },
    {
        icon: PredictiveToxicology,
        title: "Predictive Toxicology",
        description: "Deep learning and biophysical simulations predict off-target toxicity and early ocular/blood risk heatmaps with 80%+ accuracy.",
        color: "from-purple-500 to-pink-500"
    },
    {
        icon: IPRiskReasoning,
        title: "IP Risk Reasoning",
        description: "AI-driven structural novelty analysis. Logic-based reasoning to navigate the patent landscape and ensure design freedom.",
        color: "from-green-500 to-emerald-500"
    },
    {
        icon: CompetitiveIntelligence,
        title: "Competitive Intelligence",
        description: "Dynamic Relative Scoring (%) against industry gold standards like Enhertu and Trodelvy via real-time benchmarking.",
        color: "from-orange-500 to-red-500"
    },
    {
        icon: TheLibrarian,
        title: "The Librarian",
        description: "Real-time mapping of scientific evidence (PMID) across 30M+ papers to provide a solid logical foundation for every candidate.",
        color: "from-indigo-500 to-violet-500"
    },
    {
        icon: InSilico3DModeling,
        title: "In Silico 3D Modeling",
        description: "Next-Gen engine for antibody-drug binding structure simulation and affinity prediction (AlphaFold 3 Integration - Roadmap).",
        color: "from-blue-500 to-cyan-500"
    }
]

export function FeatureSection() {
    return (
        <section id="features" className="py-32 relative">
            <div className="text-center mb-20">
                <h2 className="text-3xl md:text-5xl font-bold text-white mb-6">
                    Powerful AI Agents
                </h2>
                <p className="text-gray-400 text-lg max-w-2xl mx-auto">
                    6 Core Engines Accelerating Every Stage of Research
                </p>
            </div>

            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8">
                {features.map((feature, index) => (
                    <motion.div
                        key={index}
                        initial={{ opacity: 0, y: 20 }}
                        whileInView={{ opacity: 1, y: 0 }}
                        viewport={{ once: true }}
                        transition={{ delay: index * 0.1 }}
                        className="group relative h-[320px] perspective-1000"
                    >
                        <div className="relative h-full w-full transition-all duration-500 transform-style-3d group-hover:rotate-x-2 group-hover:-translate-y-4">
                            {/* Card Background with Glassmorphism & Gradient Border */}
                            <div className="absolute inset-0 rounded-2xl bg-[#0F172A]/70 backdrop-blur-[20px] border border-white/10 shadow-2xl overflow-hidden">
                                {/* Top Gradient Bar */}
                                <div className={`h-2 w-full bg-gradient-to-r ${feature.color}`} />

                                <div className="p-8 h-full flex flex-col">
                                    {/* Icon Container - Floats on hover */}
                                    <div className={`w-16 h-16 rounded-2xl bg-gradient-to-br ${feature.color} p-0.5 mb-6 shadow-lg transform transition-transform duration-500 group-hover:scale-110 group-hover:translate-z-10 animate-float`}>
                                        <div className="w-full h-full bg-[#0F172A] rounded-[14px] flex items-center justify-center overflow-hidden">
                                            <img src={feature.icon} alt={feature.title} className="w-full h-full object-cover" />
                                        </div>
                                    </div>

                                    {/* Content */}
                                    <h3 className="text-xl font-bold text-white mb-3 group-hover:text-transparent group-hover:bg-clip-text group-hover:bg-gradient-to-r group-hover:from-white group-hover:to-gray-300 transition-colors">
                                        {feature.title}
                                    </h3>
                                    <p className="text-gray-400 leading-relaxed text-sm">
                                        {feature.description}
                                    </p>
                                </div>

                                {/* Hover Glow Effect */}
                                <div className={`absolute inset-0 bg-gradient-to-br ${feature.color} opacity-0 group-hover:opacity-5 transition-opacity duration-500 pointer-events-none`} />
                            </div>
                        </div>

                        {/* Shadow for 3D effect */}
                        <div className="absolute -bottom-4 left-4 right-4 h-4 bg-black/40 blur-xl rounded-[100%] opacity-0 group-hover:opacity-100 transition-opacity duration-500" />
                    </motion.div>
                ))}
            </div>
        </section>
    )
}
