import { Link, useLocation } from 'react-router-dom'
import {
    LayoutDashboard,
    Users,
    Database,
    Sparkles,
    Settings,
    ChevronLeft,
    ChevronRight,
    Shield,
    FlaskConical
} from 'lucide-react'
import { Button } from '@/components/ui/button'
import { cn } from '@/lib/utils'

interface AdminSidebarProps {
    isCollapsed: boolean
    setIsCollapsed: (collapsed: boolean) => void
}

const menuItems = [
    { icon: LayoutDashboard, label: 'Overview', href: '/admin' },
    { icon: Database, label: 'Data Operations', href: '/admin/data-operations' },
    { icon: FlaskConical, label: 'Design Runs', href: '/admin/design-runs' },
    { icon: Users, label: 'User Operations', href: '/admin/user-operations' },
    { icon: Sparkles, label: 'AI Tuning', href: '/admin/ai-tuning' },
    { icon: Settings, label: 'Settings', href: '/admin/settings' },
]

export function AdminSidebar({ isCollapsed, setIsCollapsed }: AdminSidebarProps) {
    const location = useLocation()

    return (
        <aside
            className={cn(
                "bg-slate-900 border-r border-slate-800 flex flex-col transition-all duration-300",
                isCollapsed ? "w-16" : "w-60"
            )}
        >
            {/* Logo */}
            <div className="h-16 flex items-center px-4 border-b border-slate-800">
                <Link to="/admin" className="flex items-center gap-2">
                    <div className="w-8 h-8 bg-gradient-to-br from-purple-500 to-pink-500 rounded-lg flex items-center justify-center flex-shrink-0">
                        <Shield className="w-4 h-4 text-white" />
                    </div>
                    {!isCollapsed && (
                        <span className="font-semibold text-white">Admin Console</span>
                    )}
                </Link>
            </div>

            {/* Menu */}
            <nav className="flex-1 py-4 px-2">
                <ul className="space-y-1">
                    {menuItems.map((item) => {
                        const isActive = location.pathname === item.href ||
                            (item.href !== '/admin' && location.pathname.startsWith(item.href))
                        return (
                            <li key={item.href}>
                                <Link
                                    to={item.href}
                                    className={cn(
                                        "flex items-center gap-3 px-3 py-2.5 rounded-lg transition-colors",
                                        isActive
                                            ? "bg-gradient-to-r from-purple-500/20 to-pink-500/20 text-white border border-purple-500/30"
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
    )
}
