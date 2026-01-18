import { useState } from 'react'
import { Link, useNavigate } from 'react-router-dom'
import { motion } from 'framer-motion'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Label } from '@/components/ui/label'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Separator } from '@/components/ui/separator'
import { Checkbox } from '@/components/ui/checkbox'
import { FlaskConical, Mail, Lock, User, Building2, Loader2, ArrowRight, Chrome, Check } from 'lucide-react'
import { toast } from 'sonner'
import { signUpWithEmail, signInWithGoogle } from '@/lib/supabase'

export function SignupPage() {
    const navigate = useNavigate()
    const [formData, setFormData] = useState({
        name: '',
        email: '',
        password: '',
        confirmPassword: '',
        organization: '',
        agreeTerms: false
    })
    const [isLoading, setIsLoading] = useState(false)

    const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        const { name, value } = e.target
        setFormData(prev => ({ ...prev, [name]: value }))
    }

    const handleSignup = async (e: React.FormEvent) => {
        e.preventDefault()

        // Validation
        if (formData.password !== formData.confirmPassword) {
            toast.error('비밀번호가 일치하지 않습니다.')
            return
        }

        if (formData.password.length < 8) {
            toast.error('비밀번호는 8자 이상이어야 합니다.')
            return
        }

        if (!formData.agreeTerms) {
            toast.error('이용약관에 동의해주세요.')
            return
        }

        setIsLoading(true)

        try {
            const { data, error } = await signUpWithEmail(
                formData.email,
                formData.password,
                formData.name
            )

            if (error) {
                toast.error(error.message || '회원가입에 실패했습니다.')
                return
            }

            if (data?.user) {
                toast.success('회원가입 완료! 이메일을 확인해주세요.')
                navigate('/login')
            }
        } catch (error) {
            toast.error('회원가입에 실패했습니다. 다시 시도해주세요.')
        } finally {
            setIsLoading(false)
        }
    }

    const handleGoogleSignup = async () => {
        const { error } = await signInWithGoogle()
        if (error) {
            toast.error(error.message || 'Google 회원가입에 실패했습니다.')
        }
    }

    // Password strength indicator
    const getPasswordStrength = () => {
        const password = formData.password
        if (password.length === 0) return { level: 0, text: '', color: '' }
        if (password.length < 6) return { level: 1, text: '약함', color: 'bg-red-500' }
        if (password.length < 8) return { level: 2, text: '보통', color: 'bg-yellow-500' }
        if (password.match(/^(?=.*[a-z])(?=.*[A-Z])(?=.*\d)/)) {
            return { level: 4, text: '강함', color: 'bg-green-500' }
        }
        return { level: 3, text: '양호', color: 'bg-blue-500' }
    }

    const passwordStrength = getPasswordStrength()

    return (
        <div className="min-h-screen bg-gradient-to-br from-slate-50 to-blue-50 flex items-center justify-center p-4">
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                className="w-full max-w-md"
            >
                {/* Logo */}
                <div className="text-center mb-8">
                    <Link to="/" className="inline-flex items-center gap-2">
                        <div className="w-10 h-10 bg-[#007AFF] rounded-xl flex items-center justify-center">
                            <FlaskConical className="w-6 h-6 text-white" />
                        </div>
                        <span className="text-2xl font-bold text-gray-900">ADC-GenAI</span>
                    </Link>
                </div>

                <Card className="shadow-xl">
                    <CardHeader className="text-center">
                        <CardTitle className="text-2xl">회원가입</CardTitle>
                        <CardDescription>
                            무료로 시작하고 50 크레딧을 받으세요
                        </CardDescription>
                    </CardHeader>
                    <CardContent className="space-y-4">
                        {/* Google Signup */}
                        <Button
                            variant="outline"
                            className="w-full"
                            onClick={handleGoogleSignup}
                        >
                            <Chrome className="w-5 h-5 mr-2" />
                            Google로 시작하기
                        </Button>

                        <div className="relative">
                            <Separator />
                            <span className="absolute left-1/2 top-1/2 -translate-x-1/2 -translate-y-1/2 bg-white px-2 text-sm text-gray-500">
                                또는
                            </span>
                        </div>

                        {/* Signup Form */}
                        <form onSubmit={handleSignup} className="space-y-4">
                            {/* Name */}
                            <div className="space-y-2">
                                <Label htmlFor="name">이름</Label>
                                <div className="relative">
                                    <User className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-gray-400" />
                                    <Input
                                        id="name"
                                        name="name"
                                        placeholder="홍길동"
                                        value={formData.name}
                                        onChange={handleChange}
                                        className="pl-10"
                                        required
                                    />
                                </div>
                            </div>

                            {/* Email */}
                            <div className="space-y-2">
                                <Label htmlFor="email">이메일</Label>
                                <div className="relative">
                                    <Mail className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-gray-400" />
                                    <Input
                                        id="email"
                                        name="email"
                                        type="email"
                                        placeholder="your@email.com"
                                        value={formData.email}
                                        onChange={handleChange}
                                        className="pl-10"
                                        required
                                    />
                                </div>
                            </div>

                            {/* Organization */}
                            <div className="space-y-2">
                                <Label htmlFor="organization">소속 (선택)</Label>
                                <div className="relative">
                                    <Building2 className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-gray-400" />
                                    <Input
                                        id="organization"
                                        name="organization"
                                        placeholder="회사/연구기관"
                                        value={formData.organization}
                                        onChange={handleChange}
                                        className="pl-10"
                                    />
                                </div>
                            </div>

                            {/* Password */}
                            <div className="space-y-2">
                                <Label htmlFor="password">비밀번호</Label>
                                <div className="relative">
                                    <Lock className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-gray-400" />
                                    <Input
                                        id="password"
                                        name="password"
                                        type="password"
                                        placeholder="8자 이상"
                                        value={formData.password}
                                        onChange={handleChange}
                                        className="pl-10"
                                        required
                                    />
                                </div>
                                {formData.password && (
                                    <div className="flex items-center gap-2">
                                        <div className="flex-1 h-1 bg-gray-200 rounded-full overflow-hidden">
                                            <div
                                                className={`h-full ${passwordStrength.color} transition-all`}
                                                style={{ width: `${passwordStrength.level * 25}%` }}
                                            />
                                        </div>
                                        <span className="text-xs text-gray-500">{passwordStrength.text}</span>
                                    </div>
                                )}
                            </div>

                            {/* Confirm Password */}
                            <div className="space-y-2">
                                <Label htmlFor="confirmPassword">비밀번호 확인</Label>
                                <div className="relative">
                                    <Lock className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-gray-400" />
                                    <Input
                                        id="confirmPassword"
                                        name="confirmPassword"
                                        type="password"
                                        placeholder="비밀번호 재입력"
                                        value={formData.confirmPassword}
                                        onChange={handleChange}
                                        className="pl-10"
                                        required
                                    />
                                    {formData.confirmPassword && formData.password === formData.confirmPassword && (
                                        <Check className="absolute right-3 top-1/2 -translate-y-1/2 w-4 h-4 text-green-500" />
                                    )}
                                </div>
                            </div>

                            {/* Terms */}
                            <div className="flex items-start gap-2">
                                <Checkbox
                                    id="agreeTerms"
                                    checked={formData.agreeTerms}
                                    onCheckedChange={(checked) =>
                                        setFormData(prev => ({ ...prev, agreeTerms: checked as boolean }))
                                    }
                                />
                                <label htmlFor="agreeTerms" className="text-sm text-gray-600 leading-tight">
                                    <Link to="/terms" className="text-[#007AFF] hover:underline">이용약관</Link> 및{' '}
                                    <Link to="/privacy" className="text-[#007AFF] hover:underline">개인정보처리방침</Link>에
                                    동의합니다.
                                </label>
                            </div>

                            <Button
                                type="submit"
                                className="w-full bg-[#007AFF] hover:bg-[#0066DD]"
                                disabled={isLoading}
                            >
                                {isLoading ? (
                                    <>
                                        <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                                        가입 처리 중...
                                    </>
                                ) : (
                                    <>
                                        회원가입
                                        <ArrowRight className="w-4 h-4 ml-2" />
                                    </>
                                )}
                            </Button>
                        </form>

                        {/* Benefits */}
                        <div className="bg-blue-50 rounded-lg p-3 space-y-1">
                            <p className="text-sm font-medium text-blue-900">🎁 가입 혜택</p>
                            <ul className="text-xs text-blue-700 space-y-1">
                                <li>• 50 무료 크레딧 (Deep 분석 5회)</li>
                                <li>• Golden Set 라이브러리 열람</li>
                                <li>• PDF 리포트 다운로드</li>
                            </ul>
                        </div>

                        <p className="text-center text-sm text-gray-500">
                            이미 계정이 있으신가요?{' '}
                            <Link to="/login" className="text-[#007AFF] hover:underline font-medium">
                                로그인
                            </Link>
                        </p>
                    </CardContent>
                </Card>
            </motion.div>
        </div>
    )
}
