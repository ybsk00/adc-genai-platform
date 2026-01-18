import { GNB } from './GNB'
import { HeroSection } from './HeroSection'
import { FeatureSection } from './FeatureSection'
import { HowItWorksSection } from './HowItWorksSection'
import { PricingSection } from './PricingSection'
import { Footer } from './Footer'

export function LandingPage() {
    return (
        <div className="min-h-screen bg-[#0F172A]">
            <GNB />
            <HeroSection />
            <FeatureSection />
            <HowItWorksSection />
            <PricingSection />
            <Footer />
        </div>
    )
}
