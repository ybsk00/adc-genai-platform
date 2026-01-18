import { motion } from 'framer-motion'
import { Button } from '@/components/ui/button'
import { Card, CardContent, CardDescription, CardFooter, CardHeader, CardTitle } from '@/components/ui/card'
import { Check } from 'lucide-react'
import { Link } from 'react-router-dom'

const plans = [
    {
        name: 'Starter',
        price: 'Free',
        period: '/ month',
        description: '개인 연구자를 위한 기본 플랜',
        features: [
            '기본 구조 분석',
            'Golden Set 열람 (제한적)',
            '월 3회 시뮬레이션',
            '커뮤니티 지원'
        ],
        cta: 'Start Research',
        popular: false,
        href: '/signup'
    },
    {
        name: 'Researcher',
        price: '$499',
        period: '/ report',
        description: '소규모 랩 및 바이오 벤처',
        features: [
            '심층 독성 예측 (Toxicity)',
            '전체 AI 에이전트 사용',
            '특허 침해 분석 (FTO)',
            '상세 PDF 리포트 생성',
            '이메일 지원'
        ],
        cta: 'Get Started',
        popular: true,
        href: '/pricing'
    },
    {
        name: 'Enterprise',
        price: 'Custom',
        period: '',
        description: '제약사 및 대규모 연구 기관',
        features: [
            'API 연동 및 통합',
            '전용 클라우드 서버 (Private)',
            'SLA 및 보안 계약',
            '전담 매니저 배정',
            'On-premise 옵션'
        ],
        cta: 'Contact Sales',
        popular: false,
        href: '/contact'
    }
]

export function PricingSection() {
    return (
        <section id="pricing" className="py-32 bg-[#0F172A]">
            <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
                {/* Header */}
                <motion.div
                    initial={{ opacity: 0, y: 20 }}
                    whileInView={{ opacity: 1, y: 0 }}
                    viewport={{ once: true }}
                    className="text-center mb-20"
                >
                    <h2 className="text-3xl md:text-5xl font-bold text-white mb-6 tracking-tight">
                        Flexible Research Plans
                    </h2>
                    <p className="text-lg md:text-xl text-gray-400 max-w-2xl mx-auto leading-relaxed tracking-wide">
                        연구 단계에 맞는 최적의 플랜을 선택하세요.
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
                            <Card className={`relative h-full bg-white/5 backdrop-blur-md border-white/10 ${plan.popular ? 'border-[#06B6D4] border-2 shadow-lg shadow-[#06B6D4]/20' : 'border'}`}>
                                {plan.popular && (
                                    <div className="absolute -top-3 left-1/2 -translate-x-1/2 px-3 py-1 bg-[#06B6D4] text-white text-sm rounded-full font-medium">
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
                                        <span className="text-gray-400 ml-1 text-sm">{plan.period}</span>
                                    </div>
                                    <ul className="space-y-4">
                                        {plan.features.map((feature) => (
                                            <li key={feature} className="flex items-start gap-3">
                                                <Check className="w-5 h-5 text-[#10B981] flex-shrink-0 mt-0.5" />
                                                <span className="text-gray-300 text-sm">{feature}</span>
                                            </li>
                                        ))}
                                    </ul>
                                </CardContent>
                                <CardFooter className="mt-auto">
                                    <Button
                                        className={`w-full ${plan.popular ? 'bg-[#06B6D4] hover:bg-[#0891B2] text-white' : 'bg-white/10 hover:bg-white/20 text-white border-white/20'}`}
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

