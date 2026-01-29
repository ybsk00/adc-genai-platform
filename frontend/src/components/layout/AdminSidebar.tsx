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
    FlaskConical,
    Beaker,
    Cpu,
    Wrench,
    ScrollText,
    BarChart3,
    Zap
} from 'lucide-react'
import { Button } from '@/components/ui/button'
import { cn } from '@/lib/utils'

interface AdminSidebarProps {
    isCollapsed: boolean
    setIsCollapsed: (collapsed: boolean) => void
}

interface MenuItem {
    icon: React.ComponentType<{ className?: string }>
    label: string
    href: string
    badge?: string
    highlight?: boolean
}

interface MenuSection {
    title?: string
    items: MenuItem[]
}

const menuSections: MenuSection[] = [
    {
        items: [
            { icon: LayoutDashboard, label: 'Overview', href: '/admin' },
            { icon: Database, label: 'Total Inventory', href: '/admin/inventory' },
            { icon: Beaker, label: 'ADC Total Inventory', href: '/admin/adc-inventory' },
            { icon: FlaskConical, label: 'Staging Area', href: '/admin/staging' },
            { icon: Sparkles, label: 'Golden Set Library', href: '/admin/goldenset' },
            { icon: Users, label: 'User Operations', href: '/admin/user-operations' },
        ]
    },
    {
        title: 'v2.2 Management',
        items: [
            { icon: Cpu, label: 'AI Engine Control', href: '/admin/engine-control', badge: 'NEW', highlight: true },
            { icon: Wrench, label: 'Data Refinement', href: '/admin/refinement-hub', badge: 'NEW', highlight: true },
            { icon: ScrollText, label: 'Audit & Lineage', href: '/admin/audit-lineage', badge: 'NEW', highlight: true },
            { icon: BarChart3, label: 'Analytics & Cost', href: '/admin/analytics', badge: 'NEW', highlight: true },
        ]
    },
    {
        items: [
            { icon: Settings, label: 'Settings', href: '/admin/settings' },
        ]
    }
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
            <nav className="flex-1 py-4 px-2 overflow-y-auto">
                {menuSections.map((section, sectionIndex) => (
                    <div key={sectionIndex} className={sectionIndex > 0 ? "mt-4" : ""}>
                        {/* Section Title */}
                        {section.title && !isCollapsed && (
                            <div className="px-3 py-2 mb-1">
                                <span className="text-xs font-semibold text-purple-400 uppercase tracking-wider flex items-center gap-1">
                                    <Zap className="w-3 h-3" />
                                    {section.title}
                                </span>
                            </div>
                        )}
                        {section.title && isCollapsed && (
                            <div className="mx-2 my-2 border-t border-purple-500/30" />
                        )}

                        <ul className="space-y-1">
                            {section.items.map((item) => {
                                const isActive = location.pathname === item.href ||
                                    (item.href !== '/admin' && location.pathname.startsWith(item.href))
                                return (
                                    <li key={item.href}>
                                        <Link
                                            to={item.href}
                                            className={cn(
                                                "flex items-center gap-3 px-3 py-2.5 rounded-lg transition-colors relative",
                                                isActive
                                                    ? "bg-gradient-to-r from-purple-500/20 to-pink-500/20 text-white border border-purple-500/30"
                                                    : item.highlight
                                                    ? "text-purple-300 hover:text-white hover:bg-purple-900/30 bg-purple-900/10"
                                                    : "text-slate-400 hover:text-white hover:bg-slate-800"
                                            )}
                                        >
                                            <item.icon className={cn(
                                                "w-5 h-5 flex-shrink-0",
                                                item.highlight && !isActive && "text-purple-400"
                                            )} />
                                            {!isCollapsed && (
                                                <span className="flex-1">{item.label}</span>
                                            )}
                                            {!isCollapsed && item.badge && (
                                                <span className="px-1.5 py-0.5 text-[10px] font-bold bg-gradient-to-r from-purple-600 to-pink-600 text-white rounded">
                                                    {item.badge}
                                                </span>
                                            )}
                                        </Link>
                                    </li>
                                )
                            })}
                        </ul>
                    </div>
                ))}
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
