import { useState, useEffect } from 'react'
import { Outlet, useNavigate } from 'react-router-dom'
import {
    LogOut,
    Shield,
    Loader2
} from 'lucide-react'
import { Button } from '@/components/ui/button'
import { Avatar, AvatarFallback } from '@/components/ui/avatar'
import { supabase } from '@/lib/supabase'
import { toast } from 'sonner'
import { AdminSidebar } from '@/components/layout/AdminSidebar'

/**
 * Admin Console Layout
 * 디자인: Dark Theme (Slate-900 배경)
 * 참고: 03_Frontend_Admin_Console.md
 */
export function AdminLayout() {
    const [isCollapsed, setIsCollapsed] = useState(false)
    const [isLoading, setIsLoading] = useState(true)
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
            {/* Sidebar Component */}
            <AdminSidebar isCollapsed={isCollapsed} setIsCollapsed={setIsCollapsed} />

            {/* Main Content */}
            <div className="flex-1 flex flex-col">
                {/* Top Header */}
                <header className="h-16 bg-slate-900 border-b border-slate-800 flex items-center justify-between px-6">
                    {/* Breadcrumb */}
                    <div className="text-sm">
                        <span className="text-slate-500">Admin</span>
                        <span className="mx-2 text-slate-600">/</span>
                        <span className="text-white font-medium">Console</span>
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
                            <Button
                                variant="ghost"
                                size="icon"
                                className="text-slate-400 hover:text-white"
                                onClick={async () => {
                                    await supabase?.auth.signOut()
                                    toast.success('로그아웃 되었습니다.')
                                    navigate('/login')
                                }}
                            >
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

