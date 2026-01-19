import { useState, useEffect } from 'react'
import { Link, Outlet, useLocation, useNavigate } from 'react-router-dom'
import {
    LayoutDashboard,
    Users,
    Database,
    Sparkles,
    Settings,
    ChevronLeft,
    ChevronRight,
    LogOut,
    Shield,
    Loader2
} from 'lucide-react'
import { Button } from '@/components/ui/button'
import { Avatar, AvatarFallback } from '@/components/ui/avatar'
import { cn } from '@/lib/utils'
import { supabase } from '@/lib/supabase'
import { toast } from 'sonner'

const menuItems = [
    { icon: LayoutDashboard, label: 'Overview', href: '/admin' },
    { icon: Users, label: 'User Operations', href: '/admin/users' },
    { icon: Database, label: 'Data Operations', href: '/admin/data' },
    { icon: Sparkles, label: 'AI Tuning', href: '/admin/ai' },
    { icon: Settings, label: 'Settings', href: '/admin/settings' },
]

/**
 * Admin Console Layout
 * 디자인: Dark Theme (Slate-900 배경)
 * 참고: 03_Frontend_Admin_Console.md
 */
export function AdminLayout() {
    const [isCollapsed, setIsCollapsed] = useState(false)
    const [isLoading, setIsLoading] = useState(true)
    const location = useLocation()
    const navigate = useNavigate()

    useEffect(() => {
        const checkAdmin = async () => {
            if (!supabase) {
                toast.error('인증 시스템이 설정되지 않았습니다.')
                navigate('/dashboard')
                return
            }

            try {
                const { data: { user } } = await supabase.auth.getUser()

                if (!user) {
                    toast.error('로그인이 필요합니다.')
                    navigate('/login')
                    return
                }

                if (user.email !== 'admin@admin.com') {
                    toast.error('관리자 권한이 없습니다.')
                    navigate('/dashboard')
                    return
                }

                setIsLoading(false)
            } catch (error) {
                console.error('Admin check error:', error)
                navigate('/dashboard')
            }
        }

        checkAdmin()
    }, [navigate])

    if (isLoading) {
        return (
            <div className="min-h-screen bg-slate-950 flex items-center justify-center">
                <Loader2 className="w-8 h-8 text-purple-500 animate-spin" />
            </div>
        )
    }

    return (
        <div className="min-h-screen bg-slate-950 flex">
            {/* Dark Sidebar */}
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

            {/* Main Content */}
            <div className="flex-1 flex flex-col">
                {/* Top Header */}
                <header className="h-16 bg-slate-900 border-b border-slate-800 flex items-center justify-between px-6">
                    {/* Breadcrumb */}
                    <div className="text-sm">
                        <span className="text-slate-500">Admin</span>
                        <span className="mx-2 text-slate-600">/</span>
                        <span className="text-white font-medium">Overview</span>
                    </div>

                    {/* Right Section */}
                    <div className="flex items-center gap-4">
                        {/* Admin Badge */}
                        <div className="flex items-center gap-2 px-3 py-1.5 bg-purple-500/20 text-purple-400 rounded-full text-sm border border-purple-500/30">
                            <Shield className="w-4 h-4" />
                            <span className="font-medium">Super Admin</span>
                        </div>

                        {/* User Menu */}
                        <div className="flex items-center gap-3 pl-4 border-l border-slate-700">
                            <Avatar className="w-8 h-8">
                                <AvatarFallback className="bg-gradient-to-br from-purple-500 to-pink-500 text-white text-sm">
                                    AD
                                </AvatarFallback>
                            </Avatar>
                            <Button variant="ghost" size="icon" className="text-slate-400 hover:text-white">
                                <LogOut className="w-4 h-4" />
                            </Button>
                        </div>
                    </div>
                </header>

                {/* Page Content */}
                <main className="flex-1 p-6 overflow-auto bg-slate-950">
                    <Outlet />
                </main>
            </div>
        </div>
    )
}
