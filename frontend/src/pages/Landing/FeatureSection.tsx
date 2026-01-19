import { motion } from 'framer-motion'
import { Beaker, Database, FileText, Globe, Layers, Zap } from 'lucide-react'

const features = [
    {
        icon: Layers,
        title: "In Silico 3D Modeling",
        description: "AlphaFold3-based antibody-drug binding structure simulation and binding affinity prediction",
        color: "from-blue-500 to-cyan-500"
    },
    {
        icon: Zap,
        title: "Predictive Toxicology",
        description: "Off-target toxicity heatmap analysis and early detection of ocular/blood toxicity risks",
        color: "from-purple-500 to-pink-500"
    },
    {
        icon: Globe,
        title: "IP Landscape Analysis",
        description: "Real-time FTO traffic light system integrated with global patent databases",
        color: "from-green-500 to-emerald-500"
    },
    {
        icon: FileText,
        title: "Clinical Strategy",
        description: "Automated Phase 1 protocol generation based on FDA guidelines and clinical success probability reports",
        color: "from-indigo-500 to-violet-500"
    },
    {
        icon: Database,
        title: "Market Intelligence",
        description: "Real-time tracking of global pipelines and competitor development status dashboard",
        color: "from-orange-500 to-red-500"
    },
    {
        icon: Beaker,
        title: "Golden Set Library",
        description: "Comparative analysis library of FDA-approved reference drug mechanism data vs. in-house experimental data",
        color: "from-yellow-500 to-amber-500"
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
                            <div className="absolute inset-0 rounded-2xl bg-[#1E293B]/80 backdrop-blur-xl border border-white/10 shadow-2xl overflow-hidden">
                                {/* Top Gradient Bar */}
                                <div className={`h-2 w-full bg-gradient-to-r ${feature.color}`} />

                                <div className="p-8 h-full flex flex-col">
                                    {/* Icon Container - Floats on hover */}
                                    <div className={`w-16 h-16 rounded-2xl bg-gradient-to-br ${feature.color} p-0.5 mb-6 shadow-lg transform transition-transform duration-500 group-hover:scale-110 group-hover:translate-z-10`}>
                                        <div className="w-full h-full bg-[#0F172A] rounded-[14px] flex items-center justify-center">
                                            <feature.icon className="w-8 h-8 text-white" />
                                        </div>
                                    </div>

                                    {/* Content */}
                                    <h3 className="text-xl font-bold text-white mb-3 group-hover:text-transparent group-hover:bg-clip-text group-hover:bg-gradient-to-r group-hover:from-white group-hover:to-gray-300 transition-colors">
                                        {feature.title}
                                    </h3>
                                    <p className="text-gray-400 leading-relaxed text-sm">
                                        {feature.description}
                                    </p>

                                    {/* Bottom Action Hint - Removed as per request */}
                                    {/* <div className="mt-auto pt-6 flex items-center text-sm font-medium text-gray-500 group-hover:text-white transition-colors">
                                        <span>View Details</span>
                                        <div className={`ml-2 w-8 h-[1px] bg-gradient-to-r ${feature.color}`} />
                                    </div> */}
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
