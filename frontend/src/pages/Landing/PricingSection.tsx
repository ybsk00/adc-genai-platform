import { motion } from 'framer-motion'
import { Button } from '@/components/ui/button'
import { Card, CardContent, CardDescription, CardFooter, CardHeader, CardTitle } from '@/components/ui/card'
import { Check } from 'lucide-react'
import { Link } from 'react-router-dom'

const plans = [
    {
        name: 'Researcher',
        price: '$0',
        period: '무료',
        description: '바이오 연구자를 위한 기본 플랜',
        features: [
            'LIV-1 기본 분석',
            'Golden Set 라이브러리 열람',
            '월 3회 분석',
            '기본 PDF 리포트'
        ],
        cta: 'Sign Up Free',
        popular: false,
        href: '/signup'
    },
    {
        name: 'Developer',
        price: '$499',
        period: '/ 리포트',
        description: 'CDMO 및 개발팀을 위한 프로 플랜',
        features: [
            '전체 독성 예측',
            '특허 침해 검토',
            'Linker 최적화 제안',
            '6-Agent 풀가동 분석',
            '우선 처리 지원'
        ],
        cta: 'Buy Credits',
        popular: true,
        href: '/pricing'
    },
    {
        name: 'Organization',
        price: 'Custom',
        period: '문의',
        description: '대규모 팀을 위한 엔터프라이즈',
        features: [
            'API 접근 권한',
            '전용 서버',
            '무제한 시트',
            '전담 지원',
            '커스텀 통합'
        ],
        cta: 'Contact Sales',
        popular: false,
        href: '/contact'
    }
]

export function PricingSection() {
    return (
        <section id="pricing" className="py-24 bg-gradient-to-b from-[#0A0F1C] to-[#1A1F2E]">
            <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
                {/* Header */}
                <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    whileInView={{ opacity: 1, y: 0 }}
                    viewport={{ once: true }}
                    className="text-center mb-16"
                >
                    <h2 className="text-3xl md:text-4xl font-bold text-white mb-4">
                        투명한 가격 정책
                    </h2>
                    <p className="text-lg text-gray-400 max-w-2xl mx-auto">
                        필요에 맞는 플랜을 선택하세요. 언제든지 업그레이드 가능합니다.
                    </p>
                </motion.div>

                {/* Pricing Cards */}
                <div className="grid md:grid-cols-3 gap-8">
                    {plans.map((plan, index) => (
                        <motion.div
                            key={plan.name}
                            initial={{ opacity: 0, y: 20 }}
                            whileInView={{ opacity: 1, y: 0 }}
                            viewport={{ once: true }}
                            transition={{ delay: index * 0.1 }}
                        >
                            <Card className={`relative h-full bg-white/5 border-white/10 ${plan.popular ? 'border-[#007AFF] border-2 shadow-lg shadow-[#007AFF]/20' : ''}`}>
                                {plan.popular && (
                                    <div className="absolute -top-3 left-1/2 -translate-x-1/2 px-3 py-1 bg-[#007AFF] text-white text-sm rounded-full">
                                        Most Popular
                                    </div>
                                )}
                                <CardHeader>
                                    <CardTitle className="text-xl text-white">{plan.name}</CardTitle>
                                    <CardDescription className="text-gray-400">{plan.description}</CardDescription>
                                </CardHeader>
                                <CardContent>
                                    <div className="mb-6">
                                        <span className="text-4xl font-bold text-white">{plan.price}</span>
                                        <span className="text-gray-400 ml-1">{plan.period}</span>
                                    </div>
                                    <ul className="space-y-3">
                                        {plan.features.map((feature) => (
                                            <li key={feature} className="flex items-center gap-2">
                                                <Check className="w-5 h-5 text-green-500 flex-shrink-0" />
                                                <span className="text-gray-300">{feature}</span>
                                            </li>
                                        ))}
                                    </ul>
                                </CardContent>
                                <CardFooter>
                                    <Button
                                        className={`w-full ${plan.popular ? 'bg-[#007AFF] hover:bg-[#0056b3]' : 'bg-white/10 hover:bg-white/20 text-white border-white/20'}`}
                                        variant={plan.popular ? 'default' : 'outline'}
                                        asChild
                                    >
                                        <Link to={plan.href}>{plan.cta}</Link>
                                    </Button>
                                </CardFooter>
                            </Card>
                        </motion.div>
                    ))}
                </div>
            </div>
        </section>
    )
}

