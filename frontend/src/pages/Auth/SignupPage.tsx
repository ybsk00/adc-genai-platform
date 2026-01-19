import { useState } from 'react'
import { Link, useNavigate } from 'react-router-dom'
import { motion } from 'framer-motion'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Label } from '@/components/ui/label'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Separator } from '@/components/ui/separator'
import { Checkbox } from '@/components/ui/checkbox'
import { Loader2, ArrowLeft, Check } from 'lucide-react'
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
        if (password.length < 6) return { level: 1, text: 'Weak', color: 'bg-red-500' }
        if (password.length < 8) return { level: 2, text: 'Fair', color: 'bg-yellow-500' }
        if (password.match(/^(?=.*[a-z])(?=.*[A-Z])(?=.*\d)/)) {
            return { level: 4, text: 'Strong', color: 'bg-green-500' }
        }
        return { level: 3, text: 'Good', color: 'bg-blue-500' }
    }

    const passwordStrength = getPasswordStrength()

    return (
        <div className="min-h-screen bg-[#0F172A] flex items-center justify-center p-4">
            <motion.div
                initial={{ opacity: 0, scale: 0.95 }}
                animate={{ opacity: 1, scale: 1 }}
                className="w-full max-w-[400px]"
            >
                <Card className="bg-[#1E293B] border-slate-800 shadow-2xl">
                    <CardHeader className="space-y-1 pb-6">
                        <div className="flex items-center gap-2 mb-4">
                            <Link to="/" className="text-slate-400 hover:text-white transition-colors">
                                <ArrowLeft className="w-5 h-5" />
                            </Link>
                        </div>
                        <CardTitle className="text-2xl font-bold text-center text-white">Sign up</CardTitle>
                        <p className="text-center text-slate-400 text-sm">
                            Start for free and get 50 credits
                        </p>
                    </CardHeader>
                    <CardContent className="space-y-6">
                        <form onSubmit={handleSignup} className="space-y-4">
                            <div className="space-y-2">
                                <Label htmlFor="name" className="text-slate-300">Name</Label>
                                <Input
                                    id="name"
                                    name="name"
                                    placeholder="John Doe"
                                    value={formData.name}
                                    onChange={handleChange}
                                    className="bg-[#0F172A] border-slate-700 text-white placeholder:text-slate-600 h-12"
                                    required
                                />
                            </div>

                            <div className="space-y-2">
                                <Label htmlFor="email" className="text-slate-300">Email</Label>
                                <Input
                                    id="email"
                                    name="email"
                                    type="email"
                                    placeholder="your@email.com"
                                    value={formData.email}
                                    onChange={handleChange}
                                    className="bg-[#0F172A] border-slate-700 text-white placeholder:text-slate-600 h-12"
                                    required
                                />
                            </div>

                            <div className="space-y-2">
                                <Label htmlFor="password" className="text-slate-300">Password</Label>
                                <Input
                                    id="password"
                                    name="password"
                                    type="password"
                                    placeholder="8+ characters"
                                    value={formData.password}
                                    onChange={handleChange}
                                    className="bg-[#0F172A] border-slate-700 text-white placeholder:text-slate-600 h-12"
                                    required
                                />
                                {formData.password && (
                                    <div className="flex items-center gap-2 mt-2">
                                        <div className="flex-1 h-1 bg-slate-700 rounded-full overflow-hidden">
                                            <div
                                                className={`h-full ${passwordStrength.color} transition-all`}
                                                style={{ width: `${passwordStrength.level * 25}%` }}
                                            />
                                        </div>
                                        <span className="text-xs text-slate-400">{passwordStrength.text}</span>
                                    </div>
                                )}
                            </div>

                            <div className="space-y-2">
                                <Label htmlFor="confirmPassword" className="text-slate-300">Confirm Password</Label>
                                <div className="relative">
                                    <Input
                                        id="confirmPassword"
                                        name="confirmPassword"
                                        type="password"
                                        placeholder="Re-enter password"
                                        value={formData.confirmPassword}
                                        onChange={handleChange}
                                        className="bg-[#0F172A] border-slate-700 text-white placeholder:text-slate-600 h-12"
                                        required
                                    />
                                    {formData.confirmPassword && formData.password === formData.confirmPassword && (
                                        <Check className="absolute right-3 top-1/2 -translate-y-1/2 w-4 h-4 text-green-500" />
                                    )}
                                </div>
                            </div>

                            <div className="flex items-start gap-2 pt-2">
                                <Checkbox
                                    id="agreeTerms"
                                    checked={formData.agreeTerms}
                                    onCheckedChange={(checked) =>
                                        setFormData(prev => ({ ...prev, agreeTerms: checked as boolean }))
                                    }
                                    className="border-slate-600 data-[state=checked]:bg-blue-600 data-[state=checked]:border-blue-600"
                                />
                                <label htmlFor="agreeTerms" className="text-sm text-slate-400 leading-tight">
                                    I agree to the <Link to="/terms" className="text-blue-400 hover:underline">Terms of Service</Link> and <Link to="/privacy" className="text-blue-400 hover:underline">Privacy Policy</Link>.
                                </label>
                            </div>

                            <Button
                                type="submit"
                                className="w-full bg-[#3B82F6] hover:bg-[#2563EB] text-white h-12 text-base font-medium"
                                disabled={isLoading}
                            >
                                {isLoading ? <Loader2 className="w-4 h-4 animate-spin" /> : 'Sign up'}
                            </Button>
                        </form>

                        <div className="text-center text-sm">
                            <span className="text-slate-500">Already have an account? </span>
                            <Link to="/login" className="text-[#3B82F6] hover:underline font-medium">
                                Log in
                            </Link>
                        </div>

                        <div className="relative">
                            <div className="absolute inset-0 flex items-center">
                                <Separator className="bg-slate-800" />
                            </div>
                            <div className="relative flex justify-center text-xs uppercase">
                                <span className="bg-[#1E293B] px-2 text-slate-500">Or sign up with</span>
                            </div>
                        </div>

                        <div className="space-y-3">
                            <Button
                                variant="outline"
                                className="w-full bg-white hover:bg-gray-50 text-slate-700 border border-slate-200 h-12 font-medium transition-all hover:shadow-md"
                                onClick={handleGoogleSignup}
                            >
                                <svg className="w-5 h-5 mr-2" viewBox="0 0 24 24">
                                    <path
                                        fill="#4285F4"
                                        d="M22.56 12.25c0-.78-.07-1.53-.2-2.25H12v4.26h5.92c-.26 1.37-1.04 2.53-2.21 3.31v2.77h3.57c2.08-1.92 3.28-4.74 3.28-8.09z"
                                    />
                                    <path
                                        fill="#34A853"
                                        d="M12 23c2.97 0 5.46-.98 7.28-2.66l-3.57-2.77c-.98.66-2.23 1.06-3.71 1.06-2.86 0-5.29-1.93-6.16-4.53H2.18v2.84C3.99 20.53 7.7 23 12 23z"
                                    />
                                    <path
                                        fill="#FBBC05"
                                        d="M5.84 14.09c-.22-.66-.35-1.36-.35-2.09s.13-1.43.35-2.09V7.07H2.18C1.43 8.55 1 10.22 1 12s.43 3.45 1.18 4.93l2.85-2.22.81-.62z"
                                    />
                                    <path
                                        fill="#EA4335"
                                        d="M12 5.38c1.62 0 3.06.56 4.21 1.64l3.15-3.15C17.45 2.09 14.97 1 12 1 7.7 1 3.99 3.47 2.18 7.07l3.66 2.84c.87-2.6 3.3-4.53 6.16-4.53z"
                                    />
                                </svg>
                                <span className="text-slate-500 mr-1">Sign up with</span>
                                <span className="font-bold">
                                    <span className="text-[#4285F4]">G</span>
                                    <span className="text-[#EA4335]">o</span>
                                    <span className="text-[#FBBC05]">o</span>
                                    <span className="text-[#4285F4]">g</span>
                                    <span className="text-[#34A853]">l</span>
                                    <span className="text-[#EA4335]">e</span>
                                </span>
                            </Button>
                        </div>
                    </CardContent>
                </Card>
            </motion.div>
        </div>
    )
}
