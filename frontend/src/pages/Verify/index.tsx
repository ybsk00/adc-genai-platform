/**
 * Report Verification Page
 * QR 코드 스캔을 통한 보고서 무결성 검증 페이지
 */
import { useState, useEffect } from 'react'
import { useParams, useNavigate } from 'react-router-dom'
import {
    CheckCircle,
    XCircle,
    Shield,
    FileText,
    Calendar,
    Hash,
    ExternalLink,
    Loader2,
    AlertCircle
} from 'lucide-react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import { API_BASE_URL } from '@/lib/api'

interface VerificationResult {
    valid: boolean
    report_id?: string
    session_id?: string
    generated_at?: string
    page_count?: number
    chain_hash?: string
    algorithm?: string
    verification_timestamp?: string
    session_info?: {
        type?: string
        target?: string
        status?: string
    }
    compliance?: {
        standard: string
        features: string[]
    }
    error?: string
    message?: string
}

export default function VerifyReport() {
    const { hash } = useParams<{ hash: string }>()
    const navigate = useNavigate()
    const [loading, setLoading] = useState(true)
    const [result, setResult] = useState<VerificationResult | null>(null)
    const [error, setError] = useState<string | null>(null)

    useEffect(() => {
        if (!hash) {
            setError('No verification hash provided')
            setLoading(false)
            return
        }

        verifyReport(hash)
    }, [hash])

    const verifyReport = async (reportHash: string) => {
        try {
            setLoading(true)
            setError(null)

            const response = await fetch(`${API_BASE_URL}/api/report/verify/${reportHash}`)
            const data = await response.json()

            setResult(data)
        } catch (err) {
            setError('Failed to verify report. Please try again.')
            console.error('Verification error:', err)
        } finally {
            setLoading(false)
        }
    }

    if (loading) {
        return (
            <div className="min-h-screen bg-gradient-to-br from-slate-900 via-slate-800 to-slate-900 flex items-center justify-center">
                <div className="text-center">
                    <Loader2 className="w-12 h-12 text-blue-500 animate-spin mx-auto" />
                    <p className="mt-4 text-slate-300">Verifying report...</p>
                </div>
            </div>
        )
    }

    if (error) {
        return (
            <div className="min-h-screen bg-gradient-to-br from-slate-900 via-slate-800 to-slate-900 flex items-center justify-center p-4">
                <Card className="max-w-md w-full bg-slate-800/50 border-red-500/30">
                    <CardContent className="pt-6 text-center">
                        <AlertCircle className="w-16 h-16 text-red-500 mx-auto" />
                        <h1 className="mt-4 text-xl font-bold text-white">Verification Error</h1>
                        <p className="mt-2 text-slate-400">{error}</p>
                        <Button
                            className="mt-6"
                            onClick={() => navigate('/')}
                        >
                            Return Home
                        </Button>
                    </CardContent>
                </Card>
            </div>
        )
    }

    return (
        <div className="min-h-screen bg-gradient-to-br from-slate-900 via-slate-800 to-slate-900 py-12 px-4">
            <div className="max-w-2xl mx-auto">
                {/* Header */}
                <div className="text-center mb-8">
                    <div className="inline-flex items-center gap-2 px-4 py-2 bg-blue-500/20 rounded-full border border-blue-500/30 mb-4">
                        <Shield className="w-5 h-5 text-blue-400" />
                        <span className="text-blue-300 font-medium">ADC-GenAI Report Verification</span>
                    </div>
                    <h1 className="text-3xl font-bold text-white">Report Integrity Check</h1>
                </div>

                {/* Result Card */}
                <Card className={`bg-slate-800/50 border-2 ${
                    result?.valid
                        ? 'border-green-500/50'
                        : 'border-red-500/50'
                }`}>
                    <CardContent className="pt-8 pb-6">
                        {/* Status Icon & Message */}
                        <div className="text-center mb-8">
                            {result?.valid ? (
                                <>
                                    <div className="w-24 h-24 mx-auto rounded-full bg-green-500/20 flex items-center justify-center mb-4">
                                        <CheckCircle className="w-16 h-16 text-green-500" />
                                    </div>
                                    <h2 className="text-2xl font-bold text-green-400">Report Verified</h2>
                                    <p className="mt-2 text-slate-400">
                                        This report is authentic and has not been tampered with.
                                    </p>
                                </>
                            ) : (
                                <>
                                    <div className="w-24 h-24 mx-auto rounded-full bg-red-500/20 flex items-center justify-center mb-4">
                                        <XCircle className="w-16 h-16 text-red-500" />
                                    </div>
                                    <h2 className="text-2xl font-bold text-red-400">Verification Failed</h2>
                                    <p className="mt-2 text-slate-400">
                                        {result?.message || 'This report could not be verified.'}
                                    </p>
                                </>
                            )}
                        </div>

                        {/* Details */}
                        {result?.valid && (
                            <div className="space-y-4">
                                {/* Report Info */}
                                <div className="p-4 bg-slate-700/50 rounded-lg">
                                    <h3 className="text-sm font-medium text-slate-400 mb-3">Report Information</h3>
                                    <div className="grid grid-cols-2 gap-4">
                                        <div className="flex items-center gap-2">
                                            <FileText className="w-4 h-4 text-slate-500" />
                                            <span className="text-sm text-slate-300">
                                                {result.report_id || result.session_id?.slice(0, 8)}
                                            </span>
                                        </div>
                                        <div className="flex items-center gap-2">
                                            <Calendar className="w-4 h-4 text-slate-500" />
                                            <span className="text-sm text-slate-300">
                                                {result.generated_at
                                                    ? new Date(result.generated_at).toLocaleString()
                                                    : 'N/A'}
                                            </span>
                                        </div>
                                    </div>
                                </div>

                                {/* Session Info */}
                                {result.session_info && (
                                    <div className="p-4 bg-slate-700/50 rounded-lg">
                                        <h3 className="text-sm font-medium text-slate-400 mb-3">Session Details</h3>
                                        <div className="space-y-2">
                                            {result.session_info.type && (
                                                <div className="flex justify-between">
                                                    <span className="text-slate-400">Type</span>
                                                    <span className="text-slate-200">{result.session_info.type}</span>
                                                </div>
                                            )}
                                            {result.session_info.target && (
                                                <div className="flex justify-between">
                                                    <span className="text-slate-400">Target</span>
                                                    <span className="text-slate-200">{result.session_info.target}</span>
                                                </div>
                                            )}
                                            {result.session_info.status && (
                                                <div className="flex justify-between">
                                                    <span className="text-slate-400">Status</span>
                                                    <Badge variant="outline" className="text-green-400 border-green-400/30">
                                                        {result.session_info.status}
                                                    </Badge>
                                                </div>
                                            )}
                                        </div>
                                    </div>
                                )}

                                {/* Chain Hash */}
                                <div className="p-4 bg-slate-700/50 rounded-lg">
                                    <h3 className="text-sm font-medium text-slate-400 mb-3 flex items-center gap-2">
                                        <Hash className="w-4 h-4" />
                                        Digital Seal
                                    </h3>
                                    <div className="font-mono text-xs text-cyan-400 bg-slate-900 p-3 rounded break-all">
                                        {result.chain_hash}
                                    </div>
                                    <div className="mt-2 flex items-center justify-between text-xs text-slate-500">
                                        <span>Algorithm: {result.algorithm}</span>
                                        <span>{result.page_count} pages</span>
                                    </div>
                                </div>

                                {/* Compliance */}
                                {result.compliance && (
                                    <div className="p-4 bg-green-500/10 border border-green-500/30 rounded-lg">
                                        <h3 className="text-sm font-medium text-green-400 mb-2">
                                            {result.compliance.standard} Compliant
                                        </h3>
                                        <div className="flex flex-wrap gap-2">
                                            {result.compliance.features.map((feature, idx) => (
                                                <Badge key={idx} variant="outline" className="text-green-300 border-green-500/30">
                                                    {feature}
                                                </Badge>
                                            ))}
                                        </div>
                                    </div>
                                )}

                                {/* Verification Timestamp */}
                                <div className="text-center text-xs text-slate-500 mt-4">
                                    Verified at: {result.verification_timestamp
                                        ? new Date(result.verification_timestamp).toLocaleString()
                                        : new Date().toLocaleString()}
                                </div>
                            </div>
                        )}

                        {/* Actions */}
                        <div className="mt-8 flex justify-center gap-4">
                            <Button
                                variant="outline"
                                onClick={() => navigate('/')}
                                className="border-slate-600 text-slate-300 hover:bg-slate-700"
                            >
                                Return Home
                            </Button>
                            {result?.valid && result.session_id && (
                                <Button
                                    onClick={() => navigate(`/dashboard`)}
                                    className="bg-blue-600 hover:bg-blue-700"
                                >
                                    <ExternalLink className="w-4 h-4 mr-2" />
                                    Open Dashboard
                                </Button>
                            )}
                        </div>
                    </CardContent>
                </Card>

                {/* Footer */}
                <div className="mt-8 text-center text-sm text-slate-500">
                    <p>This verification page confirms the authenticity of reports</p>
                    <p>generated by the ADC-GenAI Platform.</p>
                </div>
            </div>
        </div>
    )
}
