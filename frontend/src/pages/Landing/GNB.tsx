import { useState, useEffect } from 'react'
import { Link } from 'react-router-dom'
import { Button } from '@/components/ui/button'
import { Menu, X } from 'lucide-react'

export function GNB() {
    const [isScrolled, setIsScrolled] = useState(false)
    const [isMenuOpen, setIsMenuOpen] = useState(false)

    useEffect(() => {
        const handleScroll = () => {
            setIsScrolled(window.scrollY > 20)
        }
        window.addEventListener('scroll', handleScroll)
        return () => window.removeEventListener('scroll', handleScroll)
    }, [])

    return (
        <nav className={`fixed top-0 left-0 right-0 z-50 transition-all duration-300 ${isScrolled ? 'bg-[#0F172A]/90 backdrop-blur-md border-b border-white/10' : 'bg-transparent'}`}>
            <div className="w-full max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
                <div className="flex items-center justify-between h-16">
                    {/* Logo */}
                    <Link to="/" className="flex items-center gap-2">
                        <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-blue-500 to-cyan-500 flex items-center justify-center text-white font-bold">
                            A
                        </div>
                        <span className="text-xl font-bold text-white tracking-tight">AstraForge</span>
                    </Link>

                    {/* Desktop Menu */}
                    <div className="hidden md:flex items-center gap-8">
                        <a href="#features" className="text-sm font-medium text-gray-300 hover:text-white transition-colors">Features</a>
                        <a href="#pricing" className="text-sm font-medium text-gray-300 hover:text-white transition-colors">Pricing</a>
                        <Link to="/library" className="text-sm font-medium text-gray-300 hover:text-white transition-colors">Golden Set</Link>
                    </div>

                    {/* Auth Buttons */}
                    <div className="hidden md:flex items-center gap-4">
                        <Link to="/login" className="text-sm font-medium text-gray-300 hover:text-white transition-colors">
                            Log in
                        </Link>
                        <Button size="sm" className="bg-[#007AFF] hover:bg-[#0056b3] text-white rounded-lg" asChild>
                            <Link to="/signup">Sign up</Link>
                        </Button>
                    </div>

                    {/* Mobile Menu Button */}
                    <button
                        className="md:hidden p-2 text-gray-400 hover:text-white"
                        onClick={() => setIsMenuOpen(!isMenuOpen)}
                    >
                        {isMenuOpen ? <X className="w-6 h-6" /> : <Menu className="w-6 h-6" />}
                    </button>
                </div>

                {/* Mobile Menu */}
                {isMenuOpen && (
                    <div className="md:hidden py-4 border-t border-white/10 bg-[#0F172A]">
                        <div className="flex flex-col gap-4">
                            <a href="#features" className="text-sm font-medium text-gray-300 hover:text-white px-2" onClick={() => setIsMenuOpen(false)}>Features</a>
                            <a href="#pricing" className="text-sm font-medium text-gray-300 hover:text-white px-2" onClick={() => setIsMenuOpen(false)}>Pricing</a>
                            <Link to="/library" className="text-sm font-medium text-gray-300 hover:text-white px-2" onClick={() => setIsMenuOpen(false)}>Golden Set</Link>
                            <div className="h-px bg-white/10 my-2" />
                            <Link to="/login" className="text-sm font-medium text-gray-300 hover:text-white px-2" onClick={() => setIsMenuOpen(false)}>Log in</Link>
                            <Link to="/signup" className="text-sm font-medium text-[#007AFF] hover:text-[#0056b3] px-2" onClick={() => setIsMenuOpen(false)}>Sign up</Link>
                        </div>
                    </div>
                )}
            </div>
        </nav>
    )
}
