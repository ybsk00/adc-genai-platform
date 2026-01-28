import { useEffect, useRef } from 'react'
import { HeroSection } from './HeroSection'
import { ProblemSolutionSection } from './ProblemSolutionSection'
import { FeatureSection } from './FeatureSection'
import { HowItWorksSection } from './HowItWorksSection'
import { UseCaseSection } from './UseCaseSection'
import { TechStackSection } from './TechStackSection'
import { TrustedBySection } from './TrustedBySection'
import { CTASection } from './CTASection'

export function LandingPage() {
    const containerRef = useRef<HTMLDivElement>(null)

    useEffect(() => {
        const container = containerRef.current
        if (!container) return

        const handleMouseMove = (e: MouseEvent) => {
            const rect = container.getBoundingClientRect()
            const x = e.clientX - rect.left
            const y = e.clientY - rect.top
            container.style.setProperty('--mouse-x', `${x}px`)
            container.style.setProperty('--mouse-y', `${y}px`)
        }

        container.addEventListener('mousemove', handleMouseMove)
        return () => container.removeEventListener('mousemove', handleMouseMove)
    }, [])

    return (
        <div ref={containerRef} className="relative min-h-screen bg-[#0F172A] overflow-hidden group">
            {/* Mouse Follow Light */}
            <div
                className="pointer-events-none absolute inset-0 z-0 transition-opacity duration-300 opacity-0 group-hover:opacity-100"
                style={{
                    background: `radial-gradient(600px circle at var(--mouse-x, 50%) var(--mouse-y, 50%), rgba(59, 130, 246, 0.15), transparent 40%)`
                }}
            />

            <div className="relative z-10 flex flex-col gap-0">
                <HeroSection />
                <ProblemSolutionSection />
                <FeatureSection />
                <HowItWorksSection />
                <UseCaseSection />
                <TechStackSection />
                <TrustedBySection />
                <CTASection />
            </div>
        </div>
    )
}
