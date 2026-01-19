import { useState } from 'react'
import { Link, Outlet, useLocation } from 'react-router-dom'
import {
    LayoutDashboard,
    FlaskConical,
    Library,
    FileText,
    Settings,
    ChevronLeft,
    ChevronRight,
    LogOut,
    CreditCard
} from 'lucide-react'
import { Button } from '@/components/ui/button'
import { Avatar, AvatarFallback } from '@/components/ui/avatar'
import { cn } from '@/lib/utils'

const menuItems = [
    { icon: LayoutDashboard, label: 'Dashboard', href: '/dashboard' },
    { icon: FlaskConical, label: 'ADC Builder', href: '/dashboard/builder' },
    { icon: Library, label: 'Golden Set', href: '/dashboard/library' },
    { icon: FileText, label: 'My Reports', href: '/dashboard/reports' },
    { icon: Settings, label: 'Settings', href: '/dashboard/settings' },
]

export function DashboardLayout() {
    const [isCollapsed, setIsCollapsed] = useState(false)
    const location = useLocation()

    // Debug log to verify dark mode code is running
    console.log('Rendering DashboardLayout (Dark Mode Version)')

    return (
        <div className="min-h-screen bg-slate-950 flex">
            {/* Sidebar */}
            <aside
                className={cn(
                    "bg-slate-900 border-r border-slate-800 flex flex-col transition-all duration-300",
                    isCollapsed ? "w-16" : "w-60"
                )}
            >
                {/* Logo */}
                <div className="h-16 flex items-center px-4 border-b border-slate-800">
                    <Link to="/dashboard" className="flex items-center gap-2">
                        <div className="w-8 h-8 bg-gradient-to-br from-blue-500 to-cyan-500 rounded-lg flex items-center justify-center flex-shrink-0">
                            <span className="text-white font-bold text-sm">ADC</span>
                        </div>
                        {!isCollapsed && (
                            <span className="font-semibold text-white">ADC-GenAI</span>
                        )}
                    </Link>
                </div>

                {/* Menu */}
                <nav className="flex-1 py-4 px-2">
                    <ul className="space-y-1">
                        {menuItems.map((item) => {
                            const isActive = location.pathname === item.href
                            return (
                                <li key={item.href}>
                                    <Link
                                        to={item.href}
                                        className={cn(
                                            "flex items-center gap-3 px-3 py-2 rounded-lg transition-colors",
                                            isActive
                                                ? "bg-blue-500/10 text-blue-400 border border-blue-500/20"
                                                : "text-slate-400 hover:text-white hover:bg-slate-800"
                                        )}
                                    >
                                        <item.icon className="w-5 h-5 flex-shrink-0" />
                                        {!isCollapsed && <span>{item.label}</span>}
                                    </Link>
                                </li>
                            )
                        })}
                    </ul>
                </nav>

                {/* Collapse Button */}
                <div className="p-2 border-t border-slate-800">
                    <Button
                        variant="ghost"
                        size="sm"
                        className="w-full justify-center text-slate-400 hover:text-white hover:bg-slate-800"
                        onClick={() => setIsCollapsed(!isCollapsed)}
                    >
                        {isCollapsed ? <ChevronRight className="w-4 h-4" /> : <ChevronLeft className="w-4 h-4" />}
                    </Button>
                </div>
            </aside>

            {/* Main Content */}
            <div className="flex-1 flex flex-col">
                {/* Top Header */}
                <header className="h-16 bg-slate-900 border-b border-slate-800 flex items-center justify-between px-6">
                    {/* Breadcrumb */}
                    <div className="text-sm text-slate-500">
                        <span>Home</span>
                        <span className="mx-2">/</span>
                        <span className="text-white font-medium">Dashboard</span>
                    </div>

                    {/* Right Section */}
                    <div className="flex items-center gap-4">
                        {/* Credits */}
                        <div className="flex items-center gap-2 px-3 py-1.5 bg-blue-500/10 text-blue-400 border border-blue-500/20 rounded-full text-sm">
                            <CreditCard className="w-4 h-4" />
                            <span className="font-medium">450 Credits</span>
                        </div>

                        {/* Status */}
                        <div className="flex items-center gap-2 text-sm text-slate-400">
                            <span className="w-2 h-2 bg-green-500 rounded-full animate-pulse"></span>
                            <span>System Operational</span>
                        </div>

                        {/* User Menu */}
                        <div className="flex items-center gap-3 pl-4 border-l border-slate-800">
                            <Avatar className="w-8 h-8">
                                <AvatarFallback className="bg-gradient-to-br from-blue-500 to-cyan-500 text-white text-sm">DR</AvatarFallback>
                            </Avatar>
                            <Button variant="ghost" size="icon" className="text-slate-400 hover:text-white">
                                <LogOut className="w-4 h-4" />
                            </Button>
                        </div>
                    </div>
                </header>

                {/* Page Content */}
                <main className="flex-1 p-6 overflow-auto bg-slate-950 text-slate-200">
                    <Outlet />
                </main>
            </div>
        </div>
    )
}
