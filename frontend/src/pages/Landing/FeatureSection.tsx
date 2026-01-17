import { motion } from 'framer-motion'
import { Beaker, Shield, FileSearch, Zap, Database, FileText } from 'lucide-react'

const features = [
    {
        icon: Beaker,
        title: '3D 구조 분석',
        description: 'AlphaFold/BioNeMo로 항체-약물 접합 구조를 실시간 예측하고 시각화합니다.',
        color: 'from-blue-500 to-cyan-500'
    },
    {
        icon: Shield,
        title: '독성 예측',
        description: 'RAG 기반 유사 약물 분석으로 Ocular Toxicity, Neutropenia 등 위험을 사전에 감지합니다.',
        color: 'from-red-500 to-orange-500'
    },
    {
        icon: FileSearch,
        title: '특허 분석',
        description: '4대 플랫폼(엔허투, 시젠 등) 특허 침해 여부를 자동으로 판단합니다.',
        color: 'from-green-500 to-emerald-500'
    },
    {
        icon: Zap,
        title: '경쟁사 분석',
        description: '실시간 웹 검색으로 동일 타겟 연구 중인 경쟁사를 모니터링합니다.',
        color: 'from-purple-500 to-pink-500'
    },
    {
        icon: Database,
        title: 'Golden Set',
        description: 'FDA 승인 ADC 15종 + 임상 데이터 2,000건의 검증된 데이터베이스를 제공합니다.',
        color: 'from-yellow-500 to-amber-500'
    },
    {
        icon: FileText,
        title: 'PDF 리포트',
        description: '투자심사용 프리미엄 리포트를 자동 생성하여 의사결정을 지원합니다.',
        color: 'from-indigo-500 to-violet-500'
    }
]

export function FeatureSection() {
    return (
        <section id="features" className="py-24 bg-white">
            <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
                {/* Header */}
                <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    whileInView={{ opacity: 1, y: 0 }}
                    viewport={{ once: true }}
                    className="text-center mb-16"
                >
                    <h2 className="text-3xl md:text-4xl font-bold text-gray-900 mb-4">
                        성공 확률을 높이는 6가지 핵심 기능
                    </h2>
                    <p className="text-lg text-gray-600 max-w-2xl mx-auto">
                        6명의 AI 에이전트가 동시에 분석하여 연구원의 시간을 절약합니다.
                    </p>
                </motion.div>

                {/* Feature Grid */}
                <div className="grid md:grid-cols-2 lg:grid-cols-3 gap-8">
                    {features.map((feature, index) => (
                        <motion.div
                            key={feature.title}
                            initial={{ opacity: 0, y: 20 }}
                            whileInView={{ opacity: 1, y: 0 }}
                            viewport={{ once: true }}
                            transition={{ delay: index * 0.1 }}
                            className="group p-6 rounded-2xl bg-gray-50 hover:bg-white hover:shadow-xl transition-all duration-300 border border-transparent hover:border-gray-100"
                        >
                            <div className={`w-12 h-12 rounded-xl bg-gradient-to-br ${feature.color} flex items-center justify-center mb-4 group-hover:scale-110 transition-transform`}>
                                <feature.icon className="w-6 h-6 text-white" />
                            </div>
                            <h3 className="text-xl font-semibold text-gray-900 mb-2">
                                {feature.title}
                            </h3>
                            <p className="text-gray-600">
                                {feature.description}
                            </p>
                        </motion.div>
                    ))}
                </div>
            </div>
        </section>
    )
}
