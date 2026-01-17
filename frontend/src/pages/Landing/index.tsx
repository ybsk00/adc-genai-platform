import { GNB } from './GNB'
import { HeroSection } from './HeroSection'
import { FeatureSection } from './FeatureSection'
import { PricingSection } from './PricingSection'
import { Footer } from './Footer'

export function LandingPage() {
    return (
        <div className="min-h-screen">
            <GNB />
            <HeroSection />
            <FeatureSection />
            <PricingSection />
            <Footer />
        </div>
    )
}
