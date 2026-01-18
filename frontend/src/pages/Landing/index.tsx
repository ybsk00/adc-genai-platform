import { HeroSection } from './HeroSection'
import { FeatureSection } from './FeatureSection'
import { HowItWorksSection } from './HowItWorksSection'
import { PricingSection } from './PricingSection'

export function LandingPage() {
    return (
        <div className="space-y-32 pb-32">
            <HeroSection />
            <FeatureSection />
            <HowItWorksSection />
            <PricingSection />
        </div>
    )
}
