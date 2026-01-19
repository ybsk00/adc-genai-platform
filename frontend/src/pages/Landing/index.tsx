import { HeroSection } from './HeroSection'
import { ProblemSolutionSection } from './ProblemSolutionSection'
import { FeatureSection } from './FeatureSection'
import { HowItWorksSection } from './HowItWorksSection'
import { UseCaseSection } from './UseCaseSection'
import { TechStackSection } from './TechStackSection'
import { TrustedBySection } from './TrustedBySection'
import { CTASection } from './CTASection'

export function LandingPage() {
    return (
        <div className="flex flex-col gap-0">
            <HeroSection />
            <ProblemSolutionSection />
            <FeatureSection />
            <HowItWorksSection />
            <UseCaseSection />
            <TechStackSection />
            <TrustedBySection />
            <CTASection />
        </div>
    )
}
