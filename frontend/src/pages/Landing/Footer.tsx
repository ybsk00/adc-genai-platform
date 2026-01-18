import { Link } from 'react-router-dom'
import { Container } from '@/components/layout/Container'

export function Footer() {
    return (
        <footer className="bg-[#0A0F1C] text-gray-400 py-16 border-t border-white/10 flex justify-center" style={{ display: 'flex', justifyContent: 'center' }}>
            <Container>
                <div className="grid md:grid-cols-4 gap-8">
                    {/* Brand */}
                    <div className="col-span-1">
                        <Link to="/" className="flex items-center gap-2 mb-4">
                            <div className="w-8 h-8 bg-[#007AFF] rounded-lg flex items-center justify-center">
                                <span className="text-white font-bold text-sm">ADC</span>
                            </div>
                            <span className="font-semibold text-white">ADC-GenAI</span>
                        </Link>
                        <p className="text-sm">
                            AI-driven simulation for ADC discovery. From months to minutes.
                        </p>
                    </div>

                    {/* Product */}
                    <div>
                        <h4 className="text-white font-semibold mb-4">Product</h4>
                        <ul className="space-y-2">
                            <li><a href="#features" className="hover:text-white transition-colors">Features</a></li>
                            <li><a href="#pricing" className="hover:text-white transition-colors">Pricing</a></li>
                            <li><Link to="/goldenset" className="hover:text-white transition-colors">Golden Set</Link></li>
                            <li><Link to="/docs" className="hover:text-white transition-colors">Documentation</Link></li>
                        </ul>
                    </div>

                    {/* Company */}
                    <div>
                        <h4 className="text-white font-semibold mb-4">Company</h4>
                        <ul className="space-y-2">
                            <li><Link to="/about" className="hover:text-white transition-colors">About</Link></li>
                            <li><Link to="/blog" className="hover:text-white transition-colors">Blog</Link></li>
                            <li><Link to="/careers" className="hover:text-white transition-colors">Careers</Link></li>
                            <li><Link to="/contact" className="hover:text-white transition-colors">Contact</Link></li>
                        </ul>
                    </div>

                    {/* Legal */}
                    <div>
                        <h4 className="text-white font-semibold mb-4">Legal</h4>
                        <ul className="space-y-2">
                            <li><Link to="/privacy" className="hover:text-white transition-colors">Privacy Policy</Link></li>
                            <li><Link to="/terms" className="hover:text-white transition-colors">Terms of Service</Link></li>
                            <li><Link to="/security" className="hover:text-white transition-colors">Security</Link></li>
                        </ul>
                    </div>
                </div>

                <div className="border-t border-gray-800 mt-12 pt-8 flex flex-col md:flex-row justify-between items-center">
                    <p className="text-sm">
                        © 2026 ADC-GenAI Platform. All rights reserved.
                    </p>
                    <div className="flex gap-4 mt-4 md:mt-0">
                        <a href="https://linkedin.com" className="hover:text-white transition-colors">LinkedIn</a>
                        <a href="https://twitter.com" className="hover:text-white transition-colors">Twitter</a>
                        <a href="https://github.com" className="hover:text-white transition-colors">GitHub</a>
                    </div>
                </div>
            </Container>
        </footer>
    )
}
