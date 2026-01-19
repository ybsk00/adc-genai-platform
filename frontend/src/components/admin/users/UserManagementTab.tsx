import { useState, useEffect } from 'react'
import { motion } from 'framer-motion'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Textarea } from '@/components/ui/textarea'
import { Badge } from '@/components/ui/badge'
import {
    Table,
    TableBody,
    TableCell,
    TableHead,
    TableHeader,
    TableRow,
} from '@/components/ui/table'
import {
    Dialog,
    DialogContent,
    DialogDescription,
    DialogFooter,
    DialogHeader,
    DialogTitle,
} from '@/components/ui/dialog'
import {
    Search,
    UserPlus,
    MoreHorizontal,
    CreditCard,
    Mail,
    Ban,
    Loader2,
    CheckCircle
} from 'lucide-react'
import {
    DropdownMenu,
    DropdownMenuContent,
    DropdownMenuItem,
    DropdownMenuTrigger,
} from '@/components/ui/dropdown-menu'
import { toast } from 'sonner'
import { API_BASE_URL } from '@/lib/api'

// User Interface
interface User {
    id: string
    email: string
    full_name?: string
    credits: number
    created_at: string
    role?: string
}

export function UserManagementTab() {
    const [users, setUsers] = useState<User[]>([])
    const [loading, setLoading] = useState(true)
    const [searchQuery, setSearchQuery] = useState('')

    // Credit Grant Modal State
    const [creditModalOpen, setCreditModalOpen] = useState(false)
    const [selectedUser, setSelectedUser] = useState<{ id: string; name: string } | null>(null)
    const [creditAmount, setCreditAmount] = useState('')
    const [creditReason, setCreditReason] = useState('')
    const [isGranting, setIsGranting] = useState(false)
    const [grantSuccess, setGrantSuccess] = useState(false)

    // Fetch Users
    const fetchUsers = async () => {
        try {
            setLoading(true)
            const response = await fetch(`${API_BASE_URL}/api/admin/users`)
            if (!response.ok) throw new Error('Failed to fetch users')
            const data = await response.json()
            setUsers(data)
        } catch (error) {
            console.error('Error fetching users:', error)
            toast.error('사용자 목록을 불러오는데 실패했습니다.')
        } finally {
            setLoading(false)
        }
    }

    useEffect(() => {
        fetchUsers()
    }, [])

    const filteredUsers = users.filter(user =>
        (user.full_name || '').toLowerCase().includes(searchQuery.toLowerCase()) ||
        user.email.toLowerCase().includes(searchQuery.toLowerCase())
    )

    const handleOpenCreditModal = (userId: string, userName: string) => {
        setSelectedUser({ id: userId, name: userName })
        setCreditAmount('')
        setCreditReason('')
        setGrantSuccess(false)
        setCreditModalOpen(true)
    }

    const handleGrantCredits = async () => {
        if (!selectedUser || !creditAmount || !creditReason) return

        setIsGranting(true)

        try {
            const response = await fetch(`${API_BASE_URL}/api/admin/credits/grant`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    user_id: selectedUser.id,
                    amount: parseInt(creditAmount),
                    reason: creditReason
                })
            })

            if (!response.ok) throw new Error('Failed to grant credits')

            setGrantSuccess(true)
            toast.success(`${selectedUser.name}님에게 ${creditAmount} 크레딧이 지급되었습니다.`)

            // Refresh users to show updated credits
            fetchUsers()

            setTimeout(() => {
                setCreditModalOpen(false)
                setGrantSuccess(false)
            }, 2000)

        } catch (error) {
            toast.error('크레딧 지급에 실패했습니다.')
        } finally {
            setIsGranting(false)
        }
    }

    const handleSendEmail = (_userId: string) => {
        toast.info('이메일 발송 기능은 준비 중입니다.')
    }

    const handleSuspendUser = (_userId: string) => {
        toast.warning('계정 정지 확인 모달이 필요합니다.')
    }

    return (
        <div className="space-y-6">
            {/* Search & Filters */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.1 }}
            >
                <Card className="bg-slate-900 border-slate-800">
                    <CardContent className="p-4">
                        <div className="flex items-center gap-4">
                            <div className="relative flex-1">
                                <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-slate-400" />
                                <Input
                                    placeholder="Search users by name or email..."
                                    value={searchQuery}
                                    onChange={(e) => setSearchQuery(e.target.value)}
                                    className="pl-10 bg-slate-800 border-slate-700 text-white placeholder:text-slate-500"
                                />
                            </div>
                            <Button variant="outline" className="border-slate-700 text-slate-300 hover:text-white">
                                Filters
                            </Button>
                            <Button className="bg-purple-500 hover:bg-purple-600">
                                <UserPlus className="w-4 h-4 mr-2" />
                                Add User
                            </Button>
                        </div>
                    </CardContent>
                </Card>
            </motion.div>

            {/* Users Table */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.2 }}
            >
                <Card className="bg-slate-900 border-slate-800">
                    <CardHeader>
                        <CardTitle className="text-white">All Users</CardTitle>
                        <CardDescription className="text-slate-400">
                            총 {filteredUsers.length}명의 사용자
                        </CardDescription>
                    </CardHeader>
                    <CardContent>
                        <Table>
                            <TableHeader>
                                <TableRow className="border-slate-800 hover:bg-transparent">
                                    <TableHead className="text-slate-400">User</TableHead>
                                    <TableHead className="text-slate-400">Plan</TableHead>
                                    <TableHead className="text-slate-400">Credits</TableHead>
                                    <TableHead className="text-slate-400">Status</TableHead>
                                    <TableHead className="text-slate-400">Last Active</TableHead>
                                    <TableHead className="text-slate-400 text-right">Actions</TableHead>
                                </TableRow>
                            </TableHeader>
                            <TableBody>
                                {loading ? (
                                    <TableRow>
                                        <TableCell colSpan={6} className="h-24 text-center text-slate-400">
                                            <div className="flex items-center justify-center gap-2">
                                                <Loader2 className="w-4 h-4 animate-spin" />
                                                Loading users...
                                            </div>
                                        </TableCell>
                                    </TableRow>
                                ) : filteredUsers.length === 0 ? (
                                    <TableRow>
                                        <TableCell colSpan={6} className="h-24 text-center text-slate-400">
                                            No users found.
                                        </TableCell>
                                    </TableRow>
                                ) : (
                                    filteredUsers.map((user) => (
                                        <TableRow key={user.id} className="border-slate-800 hover:bg-slate-800/50">
                                            <TableCell>
                                                <div>
                                                    <p className="font-medium text-white">{user.full_name || 'No Name'}</p>
                                                    <p className="text-sm text-slate-400">{user.email}</p>
                                                </div>
                                            </TableCell>
                                            <TableCell>
                                                <Badge
                                                    variant="outline"
                                                    className="border-slate-600 text-slate-400"
                                                >
                                                    Standard
                                                </Badge>
                                            </TableCell>
                                            <TableCell>
                                                <span className={`font-medium ${user.credits > 0 ? 'text-white' : 'text-red-400'}`}>
                                                    {user.credits}
                                                </span>
                                            </TableCell>
                                            <TableCell>
                                                <Badge
                                                    variant="outline"
                                                    className="border-green-500/30 text-green-400"
                                                >
                                                    Active
                                                </Badge>
                                            </TableCell>
                                            <TableCell className="text-slate-400">
                                                {new Date(user.created_at).toLocaleDateString()}
                                            </TableCell>
                                            <TableCell className="text-right">
                                                <DropdownMenu>
                                                    <DropdownMenuTrigger asChild>
                                                        <Button variant="ghost" size="icon" className="text-slate-400 hover:text-white">
                                                            <MoreHorizontal className="w-4 h-4" />
                                                        </Button>
                                                    </DropdownMenuTrigger>
                                                    <DropdownMenuContent align="end" className="bg-slate-800 border-slate-700">
                                                        <DropdownMenuItem
                                                            onClick={() => handleOpenCreditModal(user.id, user.full_name || user.email)}
                                                            className="text-slate-200 focus:bg-slate-700 focus:text-white"
                                                        >
                                                            <CreditCard className="w-4 h-4 mr-2" />
                                                            Grant Credits
                                                        </DropdownMenuItem>
                                                        <DropdownMenuItem
                                                            onClick={() => handleSendEmail(user.id)}
                                                            className="text-slate-200 focus:bg-slate-700 focus:text-white"
                                                        >
                                                            <Mail className="w-4 h-4 mr-2" />
                                                            Send Email
                                                        </DropdownMenuItem>
                                                        <DropdownMenuItem
                                                            onClick={() => handleSuspendUser(user.id)}
                                                            className="text-red-400 focus:bg-red-500/20 focus:text-red-400"
                                                        >
                                                            <Ban className="w-4 h-4 mr-2" />
                                                            Suspend Account
                                                        </DropdownMenuItem>
                                                    </DropdownMenuContent>
                                                </DropdownMenu>
                                            </TableCell>
                                        </TableRow>
                                    ))
                                )}
                            </TableBody>
                        </Table>
                    </CardContent>
                </Card>
            </motion.div>

            {/* Credit Grant Modal */}
            <Dialog open={creditModalOpen} onOpenChange={setCreditModalOpen}>
                <DialogContent className="bg-slate-900 border-slate-800 text-white">
                    <DialogHeader>
                        <DialogTitle>크레딧 지급</DialogTitle>
                        <DialogDescription className="text-slate-400">
                            {selectedUser?.name}님에게 크레딧을 지급합니다.
                        </DialogDescription>
                    </DialogHeader>

                    {grantSuccess ? (
                        <div className="py-8 flex flex-col items-center gap-4">
                            <motion.div
                                initial={{ scale: 0 }}
                                animate={{ scale: 1 }}
                                className="w-16 h-16 rounded-full bg-green-500/20 flex items-center justify-center"
                            >
                                <CheckCircle className="w-8 h-8 text-green-500" />
                            </motion.div>
                            <p className="text-lg font-medium text-white">지급 완료!</p>
                            <p className="text-slate-400">{creditAmount} 크레딧이 지급되었습니다.</p>
                        </div>
                    ) : (
                        <>
                            <div className="space-y-4 py-4">
                                <div>
                                    <label className="block text-sm font-medium text-slate-300 mb-2">
                                        지급량 *
                                    </label>
                                    <Input
                                        type="number"
                                        placeholder="예: 500"
                                        value={creditAmount}
                                        onChange={(e) => setCreditAmount(e.target.value)}
                                        className="bg-slate-800 border-slate-700 text-white"
                                    />
                                </div>
                                <div>
                                    <label className="block text-sm font-medium text-slate-300 mb-2">
                                        지급 사유 * <span className="text-slate-500 text-xs">(매출 정산용 기록)</span>
                                    </label>
                                    <Textarea
                                        placeholder="예: 삼성바이오 영업 미팅 선물, 베타테스터 보상..."
                                        value={creditReason}
                                        onChange={(e) => setCreditReason(e.target.value)}
                                        className="bg-slate-800 border-slate-700 text-white min-h-[100px]"
                                    />
                                </div>
                                <div className="p-3 bg-amber-500/10 border border-amber-500/30 rounded-lg">
                                    <p className="text-sm text-amber-400">
                                        ⚠️ 이 기록은 DB <code className="bg-slate-800 px-1 rounded">transactions</code> 테이블에 저장되어 매출 정산 시 사용됩니다.
                                    </p>
                                </div>
                            </div>
                            <DialogFooter>
                                <Button
                                    variant="outline"
                                    onClick={() => setCreditModalOpen(false)}
                                    className="border-slate-700 text-slate-300"
                                >
                                    취소
                                </Button>
                                <Button
                                    onClick={handleGrantCredits}
                                    disabled={!creditAmount || !creditReason || isGranting}
                                    className="bg-purple-500 hover:bg-purple-600"
                                >
                                    {isGranting ? (
                                        <>
                                            <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                                            처리 중...
                                        </>
                                    ) : (
                                        <>
                                            <CreditCard className="w-4 h-4 mr-2" />
                                            지급하기
                                        </>
                                    )}
                                </Button>
                            </DialogFooter>
                        </>
                    )}
                </DialogContent>
            </Dialog>
        </div>
    )
}
