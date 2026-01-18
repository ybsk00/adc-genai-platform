import { Beaker, Shield, FileSearch, Zap, Database, FileText } from 'lucide-react'

const features = [
    {
        icon: Beaker,
        title: 'In Silico 3D Modeling',
        description: 'AlphaFold 기반 항체-약물 결합 구조 시뮬레이션 및 안정성 예측',
        color: 'from-blue-500 to-cyan-500',
        size: 'lg:col-span-2'
    },
    {
        icon: Shield,
        title: 'Predictive Toxicology',
        description: '오프타겟(Off-target) 효과 및 독성 리스크 사전 감지 (Ocular/Neutropenia)',
        color: 'from-red-500 to-orange-500',
        size: 'lg:col-span-1'
    },
    {
        icon: FileSearch,
        title: 'IP Landscape Analysis',
        description: '엔허투, 시젠 등 주요 경쟁 약물의 특허 침해 가능성(FTO) 실시간 검토',
        color: 'from-green-500 to-emerald-500',
        size: 'lg:col-span-1'
    },
    {
        icon: Zap,
        title: 'Market Intelligence',
        description: '글로벌 파이프라인 실시간 추적 및 타겟 약물 개발 현황 모니터링',
        color: 'from-purple-500 to-pink-500',
        size: 'lg:col-span-2'
    },
    {
        icon: Database,
        title: 'Golden Set Library',
        description: 'FDA 승인 약물(Reference Drug)의 기전 데이터와 당사 실험 데이터의 비교 분석',
        color: 'from-yellow-500 to-amber-500',
        size: 'lg:col-span-1'
    },
    {
        icon: FileText,
        title: 'Automated IND Report',
        description: '의사결정을 위한 임상시험계획승인(IND) 수준의 프리미엄 리포트 자동 생성',
        color: 'from-indigo-500 to-violet-500',
        size: 'lg:col-span-2' // Adjusted to fill grid nicely if needed, or keep 1
    }
]

export function FeatureSection() {
    return (
        <section id="features" className="py-20 bg-[#0F172A]">
            <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
                {/* Header */}
                <div className="text-center mb-16">
                    <h2 className="text-3xl md:text-5xl font-bold text-white mb-6 tracking-tight">
                        Precision, Acceleration, Insight
                    </h2>
                    <p className="text-lg md:text-xl text-gray-400 max-w-2xl mx-auto leading-relaxed tracking-wide">
                        6명의 AI 에이전트가 귀사의 ADC 개발을 가속화합니다.
                    </p>
                </div>

                {/* Bento Grid */}
                <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6 auto-rows-[minmax(200px,auto)]">
                    {features.map((feature) => (
                        <div
                            key={feature.title}
                            className={`group relative p-8 rounded-3xl bg-white/5 backdrop-blur-md border border-white/10 hover:border-white/20 transition-all duration-300 hover:bg-white/10 overflow-hidden ${feature.size || ''}`}
                        >
                            {/* Hover Gradient Effect */}
                            <div className={`absolute inset-0 bg-gradient-to-br ${feature.color} opacity-0 group-hover:opacity-5 transition-opacity duration-500`} />

                            <div className="relative z-10 h-full flex flex-col justify-between">
                                <div>
                                    <div className={`w-12 h-12 rounded-2xl bg-gradient-to-br ${feature.color} flex items-center justify-center mb-6 shadow-lg`}>
                                        <feature.icon className="w-6 h-6 text-white" />
                                    </div>
                                    <h3 className="text-2xl font-bold text-white mb-3">
                                        {feature.title}
                                    </h3>
                                    <p className="text-gray-400 leading-relaxed">
                                        {feature.description}
                                    </p>
                                </div>

                                {/* Decorative Arrow */}
                                <div className="mt-6 flex justify-end opacity-0 group-hover:opacity-100 transition-opacity duration-300 transform translate-x-2 group-hover:translate-x-0">
                                    <div className="w-8 h-8 rounded-full bg-white/10 flex items-center justify-center">
                                        <svg className="w-4 h-4 text-white" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                                            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M17 8l4 4m0 0l-4 4m4-4H3" />
                                        </svg>
                                    </div>
                                </div>
                            </div>
                        </div>
                    ))}
                </div>
            </div>
        </section>
    )
}

