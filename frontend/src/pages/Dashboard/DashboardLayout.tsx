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

    return (
        <div className="min-h-screen bg-gray-50 flex">
            {/* Sidebar */}
            <aside
                className={cn(
                    "bg-white border-r border-gray-200 flex flex-col transition-all duration-300",
                    isCollapsed ? "w-16" : "w-60"
                )}
            >
                {/* Logo */}
                <div className="h-16 flex items-center px-4 border-b border-gray-100">
                    <Link to="/dashboard" className="flex items-center gap-2">
                        <div className="w-8 h-8 bg-[#007AFF] rounded-lg flex items-center justify-center flex-shrink-0">
                            <span className="text-white font-bold text-sm">ADC</span>
                        </div>
                        {!isCollapsed && (
                            <span className="font-semibold text-gray-900">ADC-GenAI</span>
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
                                                ? "bg-[#007AFF] text-white"
                                                : "text-gray-600 hover:bg-gray-100"
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
                <div className="p-2 border-t border-gray-100">
                    <Button
                        variant="ghost"
                        size="sm"
                        className="w-full justify-center"
                        onClick={() => setIsCollapsed(!isCollapsed)}
                    >
                        {isCollapsed ? <ChevronRight className="w-4 h-4" /> : <ChevronLeft className="w-4 h-4" />}
                    </Button>
                </div>
            </aside>

            {/* Main Content */}
            <div className="flex-1 flex flex-col">
                {/* Top Header */}
                <header className="h-16 bg-white border-b border-gray-200 flex items-center justify-between px-6">
                    {/* Breadcrumb */}
                    <div className="text-sm text-gray-500">
                        <span>Home</span>
                        <span className="mx-2">/</span>
                        <span className="text-gray-900 font-medium">Dashboard</span>
                    </div>

                    {/* Right Section */}
                    <div className="flex items-center gap-4">
                        {/* Credits */}
                        <div className="flex items-center gap-2 px-3 py-1.5 bg-amber-50 text-amber-700 rounded-full text-sm">
                            <CreditCard className="w-4 h-4" />
                            <span className="font-medium">450 Credits</span>
                        </div>

                        {/* Status */}
                        <div className="flex items-center gap-2 text-sm text-gray-500">
                            <span className="w-2 h-2 bg-green-500 rounded-full"></span>
                            <span>System Operational</span>
                        </div>

                        {/* User Menu */}
                        <div className="flex items-center gap-3 pl-4 border-l border-gray-200">
                            <Avatar className="w-8 h-8">
                                <AvatarFallback className="bg-[#007AFF] text-white text-sm">DR</AvatarFallback>
                            </Avatar>
                            <Button variant="ghost" size="icon">
                                <LogOut className="w-4 h-4" />
                            </Button>
                        </div>
                    </div>
                </header>

                {/* Page Content */}
                <main className="flex-1 p-6 overflow-auto">
                    <Outlet />
                </main>
            </div>
        </div>
    )
}
