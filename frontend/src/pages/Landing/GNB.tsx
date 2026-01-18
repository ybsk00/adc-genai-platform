import { useState } from 'react'
import { Link } from 'react-router-dom'
import { motion } from 'framer-motion'
import { Button } from '@/components/ui/button'
import { Menu, X } from 'lucide-react'

export function GNB() {
    const [isMenuOpen, setIsMenuOpen] = useState(false)

    return (
        <motion.nav
            initial={{ y: -100 }}
            animate={{ y: 0 }}
            className="fixed top-0 left-0 right-0 z-50 bg-[#0A0F1C]/90 backdrop-blur-md border-b border-white/10"
        >
            <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
                <div className="flex items-center justify-between h-16">
                    {/* Logo */}
                    <Link to="/" className="flex items-center gap-2">
                        <div className="w-8 h-8 bg-[#007AFF] rounded-lg flex items-center justify-center">
                            <span className="text-white font-bold text-sm">ADC</span>
                        </div>
                        <span className="font-semibold text-white">ADC-GenAI</span>
                    </Link>

                    {/* Desktop Menu */}
                    <div className="hidden md:flex items-center gap-8">
                        <a href="#features" className="text-gray-300 hover:text-white transition-colors">
                            Features
                        </a>
                        <a href="#pricing" className="text-gray-300 hover:text-white transition-colors">
                            Pricing
                        </a>
                        <a href="#goldenset" className="text-gray-300 hover:text-white transition-colors">
                            Golden Set
                        </a>
                    </div>

                    {/* Action Buttons */}
                    <div className="hidden md:flex items-center gap-3">
                        <Button variant="ghost" className="text-gray-300 hover:text-white hover:bg-white/10" asChild>
                            <Link to="/login">Log in</Link>
                        </Button>
                        <Button className="bg-[#007AFF] hover:bg-[#0056b3]" asChild>
                            <Link to="/signup">Start for Free</Link>
                        </Button>
                    </div>

                    {/* Mobile Menu Button */}
                    <button
                        className="md:hidden p-2 text-white"
                        onClick={() => setIsMenuOpen(!isMenuOpen)}
                    >
                        {isMenuOpen ? <X size={24} /> : <Menu size={24} />}
                    </button>
                </div>

                {/* Mobile Menu */}
                {isMenuOpen && (
                    <motion.div
                        initial={{ opacity: 0, height: 0 }}
                        animate={{ opacity: 1, height: 'auto' }}
                        className="md:hidden py-4 border-t border-white/10"
                    >
                        <div className="flex flex-col gap-4">
                            <a href="#features" className="text-gray-300 hover:text-white px-2 py-1">Features</a>
                            <a href="#pricing" className="text-gray-300 hover:text-white px-2 py-1">Pricing</a>
                            <a href="#goldenset" className="text-gray-300 hover:text-white px-2 py-1">Golden Set</a>
                            <div className="flex flex-col gap-2 pt-4 border-t border-white/10">
                                <Button variant="ghost" className="text-gray-300 hover:text-white hover:bg-white/10" asChild>
                                    <Link to="/login">Log in</Link>
                                </Button>
                                <Button className="bg-[#007AFF] hover:bg-[#0056b3]" asChild>
                                    <Link to="/signup">Start for Free</Link>
                                </Button>
                            </div>
                        </div>
                    </motion.div>
                )}
            </div>
        </motion.nav>
    )
}

