import { motion } from 'framer-motion'
import { ArrowRight } from 'lucide-react'

const steps = [
    {
        id: '01',
        title: 'Sequence Input',
        description: 'Input Target Antibody Sequence',
        icon: "/icons/Step_01.png",
        color: 'text-blue-400'
    },
    {
        id: '02',
        title: 'Code-Verified Reasoning',
        description: 'Simultaneous Analysis by 6 AI Agents',
        icon: "/icons/Step_02.png",
        color: 'text-purple-400'
    },
    {
        id: '03',
        title: 'Relative Benchmarking',
        description: 'Optimal Linker/Payload Recommendation',
        icon: "/icons/Step_03.png",
        color: 'text-cyan-400'
    },
    {
        id: '04',
        title: 'Certified 10-Page Report',
        description: 'Final Candidate Derivation',
        icon: "/icons/Step_04.png",
        color: 'text-green-400'
    }
]

export function HowItWorksSection() {
    return (
        <section className="relative py-24">
            <div className="text-center mb-24">
                <h2 className="text-3xl md:text-5xl font-bold text-white mb-6 tracking-tight">
                    How AstraForge Works
                </h2>
                <p className="text-lg text-gray-400 max-w-2xl mx-auto leading-relaxed">
                    From data input to results, complete in just 4 steps.
                </p>
            </div>

            <div className="relative">
                {/* Connecting Line (Desktop) */}
                <div className="hidden md:block absolute top-12 left-0 w-full h-0.5 bg-white/10">
                    <div className="absolute top-0 left-0 h-full bg-gradient-to-r from-blue-500 via-purple-500 to-green-500 w-full opacity-30" />
                </div>

                <div className="grid grid-cols-1 md:grid-cols-4 gap-8">
                    {steps.map((step, index) => (
                        <motion.div
                            key={step.id}
                            initial={{ opacity: 0, y: 20 }}
                            whileInView={{ opacity: 1, y: 0 }}
                            viewport={{ once: true }}
                            transition={{ delay: index * 0.2 }}
                            className="relative flex flex-col items-center text-center group"
                        >
                            {/* Step Number Bubble */}
                            <div className={`w-24 h-24 rounded-full bg-[#0F172A] border-2 border-white/10 group-hover:border-white/30 flex items-center justify-center mb-6 relative z-10 transition-all duration-300 group-hover:scale-110 shadow-lg shadow-black/50 overflow-hidden p-4`}>
                                <img src={step.icon} alt={step.title} className="w-full h-full object-contain" />
                                <div className="absolute -top-2 -right-2 w-8 h-8 rounded-full bg-white/10 flex items-center justify-center text-xs font-bold text-white border border-white/10">
                                    {step.id}
                                </div>
                            </div>

                            {/* Content */}
                            <h3 className="text-xl font-bold text-white mb-2">
                                {step.title}
                            </h3>
                            <p className="text-sm text-gray-400">
                                {step.description}
                            </p>

                            {/* Mobile Arrow */}
                            {index < steps.length - 1 && (
                                <div className="md:hidden mt-8 mb-4">
                                    <ArrowRight className="w-6 h-6 text-white/20 rotate-90" />
                                </div>
                            )}
                        </motion.div>
                    ))}
                </div>
            </div>
        </section>
    )
}
