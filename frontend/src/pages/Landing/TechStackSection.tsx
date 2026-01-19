import { Zap } from 'lucide-react'

export function TechStackSection() {
    return (
        <section className="py-24 bg-[#0F172A] relative overflow-hidden">
            <div className="text-center mb-16">
                <h2 className="text-3xl md:text-5xl font-bold text-white mb-6">
                    Powered by State-of-the-Art AI
                </h2>
                <p className="text-gray-400 text-lg max-w-2xl mx-auto">
                    Protect your data safely with the latest AI technology and robust security infrastructure.
                </p>
            </div>

            {/* Tech Logos */}
            <div className="flex flex-wrap justify-center gap-8 md:gap-16 mb-20 opacity-70 grayscale hover:grayscale-0 transition-all duration-500">
                <div className="flex items-center gap-2 text-2xl font-bold text-white">
                    <span className="text-[#76B900]">NVIDIA</span> BioNeMo
                </div>
                <div className="flex items-center gap-2 text-2xl font-bold text-white">
                    <span className="text-[#4285F4]">Google</span> Cloud
                </div>
                <div className="flex items-center gap-2 text-2xl font-bold text-white">
                    <span className="text-[#00C7F4]">AlphaFold</span> 3
                </div>
                <div className="flex items-center gap-2 text-2xl font-bold text-white">
                    🦜🔗 LangChain
                </div>
            </div>

            {/* Security Badges */}
            <div className="grid grid-cols-1 md:grid-cols-2 gap-6 max-w-4xl mx-auto px-4">
                <div className="flex items-center gap-4 p-6 rounded-2xl bg-slate-900 border border-white/10">
                    <div className="w-12 h-12 rounded-full bg-green-500/10 flex items-center justify-center shrink-0">
                        <Zap className="w-6 h-6 text-green-500" />
                    </div>
                    <div>
                        <h3 className="text-lg font-bold text-white">SOC2 Compliance Ready</h3>
                        <p className="text-gray-400 text-sm">Global Standard Data Security Compliance</p>
                    </div>
                </div>
                <div className="flex items-center gap-4 p-6 rounded-2xl bg-slate-900 border border-white/10">
                    <div className="w-12 h-12 rounded-full bg-blue-500/10 flex items-center justify-center shrink-0">
                        <Zap className="w-6 h-6 text-blue-500" />
                    </div>
                    <div>
                        <h3 className="text-lg font-bold text-white">On-Premise Option</h3>
                        <p className="text-gray-400 text-sm">On-premise Server Support (Data Leakage Prevention)</p>
                    </div>
                </div>
            </div>
        </section>
    )
}
