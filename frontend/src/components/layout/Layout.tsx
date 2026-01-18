import { Outlet } from 'react-router-dom'
import { GNB } from '@/pages/Landing/GNB'
import { Footer } from '@/pages/Landing/Footer'

export function Layout() {
    return (
        <div className="min-h-screen bg-[#0F172A] text-white font-sans">
            <GNB />
            <main className="w-full max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 pt-16">
                <Outlet />
            </main>
            <Footer />
        </div>
    )
}
