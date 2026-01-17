import { useState } from 'react'
import { Link, useNavigate } from 'react-router-dom'
import { motion } from 'framer-motion'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Label } from '@/components/ui/label'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Separator } from '@/components/ui/separator'
import { FlaskConical, Mail, Lock, Loader2, ArrowRight, Chrome } from 'lucide-react'
import { toast } from 'sonner'

export function LoginPage() {
    const navigate = useNavigate()
    const [email, setEmail] = useState('')
    const [password, setPassword] = useState('')
    const [isLoading, setIsLoading] = useState(false)

    const handleLogin = async (e: React.FormEvent) => {
        e.preventDefault()
        setIsLoading(true)

        try {
            // TODO: Supabase Auth 연동
            // const { data, error } = await supabase.auth.signInWithPassword({
            //   email,
            //   password,
            // })

            // Mock login
            await new Promise(resolve => setTimeout(resolve, 1000))

            toast.success('로그인 성공!')
            navigate('/dashboard')

        } catch (error) {
            toast.error('로그인에 실패했습니다. 이메일과 비밀번호를 확인해주세요.')
        } finally {
            setIsLoading(false)
        }
    }

    const handleGoogleLogin = async () => {
        // TODO: Supabase Google OAuth
        // const { data, error } = await supabase.auth.signInWithOAuth({
        //   provider: 'google',
        // })
        toast.info('Google 로그인은 준비 중입니다.')
    }

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
                        <CardTitle className="text-2xl">로그인</CardTitle>
                        <CardDescription>
                            계정에 로그인하여 ADC 분석을 시작하세요
                        </CardDescription>
                    </CardHeader>
                    <CardContent className="space-y-4">
                        {/* Google Login */}
                        <Button
                            variant="outline"
                            className="w-full"
                            onClick={handleGoogleLogin}
                        >
                            <Chrome className="w-5 h-5 mr-2" />
                            Google로 계속하기
                        </Button>

                        <div className="relative">
                            <Separator />
                            <span className="absolute left-1/2 top-1/2 -translate-x-1/2 -translate-y-1/2 bg-white px-2 text-sm text-gray-500">
                                또는
                            </span>
                        </div>

                        {/* Email Login Form */}
                        <form onSubmit={handleLogin} className="space-y-4">
                            <div className="space-y-2">
                                <Label htmlFor="email">이메일</Label>
                                <div className="relative">
                                    <Mail className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-gray-400" />
                                    <Input
                                        id="email"
                                        type="email"
                                        placeholder="your@email.com"
                                        value={email}
                                        onChange={(e) => setEmail(e.target.value)}
                                        className="pl-10"
                                        required
                                    />
                                </div>
                            </div>

                            <div className="space-y-2">
                                <div className="flex items-center justify-between">
                                    <Label htmlFor="password">비밀번호</Label>
                                    <Link
                                        to="/forgot-password"
                                        className="text-sm text-[#007AFF] hover:underline"
                                    >
                                        비밀번호 찾기
                                    </Link>
                                </div>
                                <div className="relative">
                                    <Lock className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-gray-400" />
                                    <Input
                                        id="password"
                                        type="password"
                                        placeholder="••••••••"
                                        value={password}
                                        onChange={(e) => setPassword(e.target.value)}
                                        className="pl-10"
                                        required
                                    />
                                </div>
                            </div>

                            <Button
                                type="submit"
                                className="w-full bg-[#007AFF] hover:bg-[#0066DD]"
                                disabled={isLoading}
                            >
                                {isLoading ? (
                                    <>
                                        <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                                        로그인 중...
                                    </>
                                ) : (
                                    <>
                                        로그인
                                        <ArrowRight className="w-4 h-4 ml-2" />
                                    </>
                                )}
                            </Button>
                        </form>

                        <p className="text-center text-sm text-gray-500">
                            계정이 없으신가요?{' '}
                            <Link to="/signup" className="text-[#007AFF] hover:underline font-medium">
                                회원가입
                            </Link>
                        </p>
                    </CardContent>
                </Card>

                <p className="text-center text-xs text-gray-400 mt-6">
                    로그인하면 <Link to="/terms" className="underline">이용약관</Link> 및{' '}
                    <Link to="/privacy" className="underline">개인정보처리방침</Link>에 동의하는 것으로 간주됩니다.
                </p>
            </motion.div>
        </div>
    )
}
