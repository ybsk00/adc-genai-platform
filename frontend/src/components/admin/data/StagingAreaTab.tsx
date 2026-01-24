import { useState, useEffect } from 'react'
import { motion } from 'framer-motion'
import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query'
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import { Input } from '@/components/ui/input'
import {
    Select,
    SelectContent,
    SelectItem,
    SelectTrigger,
    SelectValue,
} from '@/components/ui/select'
import {
    Dialog,
    DialogContent,
    DialogDescription,
    DialogFooter,
    DialogHeader,
    DialogTitle,
} from '@/components/ui/dialog'
import { Textarea } from '@/components/ui/textarea'
import { ScrollArea } from '@/components/ui/scroll-area'
import {
    CheckCircle, XCircle, FileText, Loader2, ArrowRight, Sparkles,
    AlertCircle, Clock, CheckCircle2, ChevronLeft, ChevronRight, Search, Brain
} from 'lucide-react'
import { toast } from 'sonner'
import { getSession } from '@/lib/supabase'
import { API_BASE_URL } from '@/lib/api'

interface GoldenSetDraft {
    id: string
    drug_name: string
    target: string
    payload?: string
    linker?: string
    smiles_code?: string
    enrichment_source: string
    created_at: string
    outcome_type?: string
    failure_reason?: string
    is_ai_extracted?: boolean
    raw_data?: string
    properties?: Record<string, unknown>
}

interface DraftsResponse {
    data: GoldenSetDraft[]
    total: number
    limit: number
    offset: number
}

export function StagingAreaTab() {
    const queryClient = useQueryClient()
    const [selectedDraft, setSelectedDraft] = useState<GoldenSetDraft | null>(null)
    const [rejectReason, setRejectReason] = useState('')
    const [isRejectDialogOpen, setIsRejectDialogOpen] = useState(false)

    // Pagination & Search State
    const [page, setPage] = useState(1)
    const [pageSize] = useState(20)
    const [searchQuery, setSearchQuery] = useState('')
    const [sourceFilter, setSourceFilter] = useState('')
    const [debouncedSearch, setDebouncedSearch] = useState('')

    // Debounce search
    useEffect(() => {
        const timer = setTimeout(() => setDebouncedSearch(searchQuery), 300)
        return () => clearTimeout(timer)
    }, [searchQuery])

    // Fetch Drafts with Pagination
    const { data: draftsResponse, isLoading, error } = useQuery({
        queryKey: ['goldenSetDrafts', page, pageSize, debouncedSearch, sourceFilter],
        queryFn: async () => {
            const { session } = await getSession()
            // 세션이 없어도 일단 진행 (API가 401 줄 것임), 디버깅 위해
            
            const offset = (page - 1) * pageSize
            const params = new URLSearchParams({
                limit: String(pageSize),
                offset: String(offset),
                ...(debouncedSearch && { search: debouncedSearch }),
                ...(sourceFilter && { source: sourceFilter })
            })

            const response = await fetch(`${API_BASE_URL}/api/admin/goldenset/drafts?${params}`, {
                headers: {
                    Authorization: `Bearer ${session?.access_token || ''}`
                }
            })
            if (!response.ok) {
                const errText = await response.text()
                throw new Error(`Failed to fetch drafts: ${response.status} ${errText}`)
            }
            return response.json() as Promise<DraftsResponse>
        },
        retry: 1
    })

    const drafts = draftsResponse?.data || []
    const totalCount = draftsResponse?.total || 0
    const totalPages = Math.ceil(totalCount / pageSize)

    // ... (mutations remain same)

    if (error) {
        return (
            <div className="flex flex-col items-center justify-center h-full text-red-400 p-8">
                <AlertCircle className="w-12 h-12 mb-4" />
                <p className="text-lg font-semibold">데이터를 불러오지 못했습니다.</p>
                <p className="text-sm text-slate-500 mt-2">{error.message}</p>
                <Button variant="outline" className="mt-4" onClick={() => queryClient.invalidateQueries({ queryKey: ['goldenSetDrafts'] })}>
                    <RefreshCw className="w-4 h-4 mr-2" /> 재시도
                </Button>
            </div>
        )
    }

    return (
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6 h-[calc(100vh-200px)]">
        mutationFn: async (id: string) => {
            const { session } = await getSession()
            if (!session) throw new Error('No session')

            const response = await fetch(`${API_BASE_URL}/api/admin/goldenset/${id}/refine`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    Authorization: `Bearer ${session.access_token}`
                }
            })
            if (!response.ok) throw new Error('Failed to refine')
            return response.json()
        },
        onSuccess: (data) => {
            if (data.status === 'success') {
                toast.success(`AI 분석 완료: Target=${data.analysis?.target || 'Unknown'}`)
                queryClient.invalidateQueries({ queryKey: ['goldenSetDrafts'] })
                // 선택된 드래프트 업데이트
                if (selectedDraft) {
                    setSelectedDraft({
                        ...selectedDraft,
                        target: data.analysis?.target || selectedDraft.target,
                        smiles_code: data.smiles_code || selectedDraft.smiles_code
                    })
                }
            } else {
                toast.error(data.message || 'AI 분석 실패')
            }
        },
        onError: () => {
            toast.error('AI 분석 요청 실패')
        }
    })

    // Approve Mutation
    const approveMutation = useMutation({
        mutationFn: async (id: string) => {
            const { session } = await getSession()
            if (!session) throw new Error('No session')

            const response = await fetch(`${API_BASE_URL}/api/admin/goldenset/${id}/approve`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    Authorization: `Bearer ${session.access_token}`
                },
                body: JSON.stringify({})
            })
            if (!response.ok) throw new Error('Failed to approve')
            return response.json()
        },
        onSuccess: () => {
            toast.success('승인되었습니다.')
            queryClient.invalidateQueries({ queryKey: ['goldenSetDrafts'] })
            setSelectedDraft(null)
        },
        onError: () => {
            toast.error('승인에 실패했습니다.')
        }
    })

    // Reject Mutation
    const rejectMutation = useMutation({
        mutationFn: async ({ id, reason }: { id: string, reason: string }) => {
            const { session } = await getSession()
            if (!session) throw new Error('No session')

            const response = await fetch(`${API_BASE_URL}/api/admin/goldenset/${id}/reject`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    Authorization: `Bearer ${session.access_token}`
                },
                body: JSON.stringify({ reason })
            })
            if (!response.ok) throw new Error('Failed to reject')
            return response.json()
        },
        onSuccess: () => {
            toast.success('반려되었습니다.')
            queryClient.invalidateQueries({ queryKey: ['goldenSetDrafts'] })
            setIsRejectDialogOpen(false)
            setRejectReason('')
            setSelectedDraft(null)
        },
        onError: () => {
            toast.error('반려에 실패했습니다.')
        }
    })

    const handleRejectClick = (draft: GoldenSetDraft) => {
        setSelectedDraft(draft)
        setIsRejectDialogOpen(true)
    }

    const handleApproveClick = (draft: GoldenSetDraft) => {
        if (confirm(`${draft.drug_name}을(를) 승인하시겠습니까?`)) {
            approveMutation.mutate(draft.id)
        }
    }

    return (
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6 h-[calc(100vh-200px)]">
            {/* Inbox List */}
            <Card className="bg-slate-900 border-slate-800 lg:col-span-1 flex flex-col">
                <CardHeader className="pb-2">
                    <CardTitle className="text-white flex items-center justify-between">
                        <span>Inbox</span>
                        <Badge variant="secondary" className="bg-purple-500/20 text-purple-300">
                            {totalCount} Total
                        </Badge>
                    </CardTitle>
                    <CardDescription className="text-slate-400">
                        AI가 추출한 검토 대기 항목
                    </CardDescription>

                    {/* Search & Filter Bar */}
                    <div className="flex gap-2 mt-3">
                        <div className="relative flex-1">
                            <Search className="absolute left-2 top-2.5 h-4 w-4 text-slate-500" />
                            <Input
                                placeholder="약물명 검색..."
                                value={searchQuery}
                                onChange={(e) => {
                                    setSearchQuery(e.target.value)
                                    setPage(1)
                                }}
                                className="pl-8 bg-slate-800 border-slate-700 text-white text-sm"
                            />
                        </div>
                        <Select value={sourceFilter} onValueChange={(v) => {
                            setSourceFilter(v === 'all' ? '' : v)
                            setPage(1)
                        }}>
                            <SelectTrigger className="w-[130px] bg-slate-800 border-slate-700 text-white text-sm">
                                <SelectValue placeholder="Source" />
                            </SelectTrigger>
                            <SelectContent className="bg-slate-900 border-slate-700">
                                <SelectItem value="all">All Sources</SelectItem>
                                <SelectItem value="open_fda_api">OpenFDA</SelectItem>
                                <SelectItem value="clinical_trials_api_v2">ClinicalTrials</SelectItem>
                            </SelectContent>
                        </Select>
                    </div>
                </CardHeader>
                <CardContent className="flex-1 overflow-hidden p-0 flex flex-col">
                    <ScrollArea className="h-full">
                        {isLoading ? (
                            <div className="flex justify-center p-8">
                                <Loader2 className="w-6 h-6 text-purple-500 animate-spin" />
                            </div>
                        ) : drafts?.length === 0 ? (
                            <div className="text-center p-8 text-slate-500">
                                대기 중인 항목이 없습니다.
                            </div>
                        ) : (
                            <div className="divide-y divide-slate-800">
                                {drafts?.map((draft) => (
                                    <div
                                        key={draft.id}
                                        className={`p-4 cursor-pointer hover:bg-slate-800/50 transition-colors ${selectedDraft?.id === draft.id ? 'bg-slate-800 border-l-2 border-purple-500' : ''
                                            }`}
                                        onClick={() => setSelectedDraft(draft)}
                                    >
                                        <div className="flex justify-between items-start mb-1">
                                            <h4 className="font-medium text-white">{draft.drug_name}</h4>
                                            <span className="text-xs text-slate-500">
                                                {new Date(draft.created_at).toLocaleDateString()}
                                            </span>
                                        </div>
                                        <p className="text-sm text-slate-400 mb-2">{draft.target}</p>
                                        <div className="flex gap-2 mt-2">
                                            <Badge variant="outline" className="text-xs border-slate-700 text-slate-500">
                                                {draft.enrichment_source}
                                            </Badge>
                                            {draft.is_ai_extracted && (
                                                <Badge variant="secondary" className="bg-purple-500/10 text-purple-400 border-purple-500/20 text-[10px] gap-1">
                                                    <Sparkles className="w-2.5 h-2.5" /> ✨ AI Extracted
                                                </Badge>
                                            )}
                                        </div>
                                    </div>
                                ))}
                            </div>
                        )}
                    </ScrollArea>

                    {/* Pagination Controls */}
                    {totalPages > 1 && (
                        <div className="flex items-center justify-between p-3 border-t border-slate-800 bg-slate-900/50">
                            <span className="text-xs text-slate-500">
                                {page} / {totalPages} 페이지
                            </span>
                            <div className="flex gap-1">
                                <Button
                                    size="sm"
                                    variant="ghost"
                                    onClick={() => setPage(p => Math.max(1, p - 1))}
                                    disabled={page === 1}
                                    className="h-7 w-7 p-0"
                                >
                                    <ChevronLeft className="h-4 w-4" />
                                </Button>
                                <Button
                                    size="sm"
                                    variant="ghost"
                                    onClick={() => setPage(p => Math.min(totalPages, p + 1))}
                                    disabled={page === totalPages}
                                    className="h-7 w-7 p-0"
                                >
                                    <ChevronRight className="h-4 w-4" />
                                </Button>
                            </div>
                        </div>
                    )}
                </CardContent>
            </Card>

            {/* Comparison View */}
            <Card className="bg-slate-900 border-slate-800 lg:col-span-2 flex flex-col">
                {selectedDraft ? (
                    <>
                        <CardHeader className="border-b border-slate-800 pb-4">
                            <div className="flex justify-between items-center">
                                <div>
                                    <CardTitle className="text-white flex items-center gap-2">
                                        {selectedDraft.drug_name}
                                        {selectedDraft.outcome_type === 'Success' && (
                                            <Badge variant="outline" className="bg-green-500/10 text-green-400 border-green-500/20 text-xs gap-1">
                                                <CheckCircle2 className="w-3 h-3" /> Success
                                            </Badge>
                                        )}
                                        {selectedDraft.outcome_type === 'Failure' && (
                                            <Badge variant="outline" className="bg-red-500/10 text-red-400 border-red-500/20 text-xs gap-1">
                                                <AlertCircle className="w-3 h-3" /> Failure
                                            </Badge>
                                        )}
                                        {selectedDraft.outcome_type === 'Ongoing' && (
                                            <Badge variant="outline" className="bg-blue-500/10 text-blue-400 border-blue-500/20 text-xs gap-1">
                                                <Clock className="w-3 h-3" /> Ongoing
                                            </Badge>
                                        )}
                                    </CardTitle>
                                    <CardDescription className="text-slate-400 flex items-center gap-2 mt-1">
                                        Source: {selectedDraft.enrichment_source}
                                        <ArrowRight className="w-3 h-3" />
                                        Target: {selectedDraft.target}
                                    </CardDescription>
                                    {selectedDraft.failure_reason && (
                                        <div className="mt-2 p-2 bg-red-500/10 border border-red-500/20 rounded text-xs text-red-300 flex items-start gap-2">
                                            <AlertCircle className="w-3 h-3 mt-0.5 shrink-0" />
                                            <span>Failure Reason: {selectedDraft.failure_reason}</span>
                                        </div>
                                    )}
                                </div>
                                <div className="flex gap-2">
                                    <Button
                                        variant="outline"
                                        className="border-red-500/30 text-red-400 hover:bg-red-500/10 hover:text-red-300"
                                        onClick={() => handleRejectClick(selectedDraft)}
                                    >
                                        <XCircle className="w-4 h-4 mr-2" />
                                        Reject
                                    </Button>
                                    <Button
                                        variant="outline"
                                        className="border-purple-500/30 text-purple-300 hover:bg-purple-500/10"
                                        onClick={() => refineMutation.mutate(selectedDraft.id)}
                                        disabled={refineMutation.isPending}
                                    >
                                        {refineMutation.isPending ? (
                                            <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                                        ) : (
                                            <Brain className="w-4 h-4 mr-2" />
                                        )}
                                        AI 분석
                                    </Button>
                                    <Button
                                        className="bg-green-500 hover:bg-green-600 text-white"
                                        onClick={() => handleApproveClick(selectedDraft)}
                                    >
                                        <CheckCircle className="w-4 h-4 mr-2" />
                                        Approve
                                    </Button>
                                </div>
                            </div>
                        </CardHeader>
                        <CardContent className="flex-1 overflow-hidden p-0">
                            <div className="grid grid-cols-2 h-full divide-x divide-slate-800">
                                {/* Raw Data View */}
                                <div className="p-4 overflow-auto bg-slate-950/50">
                                    <h4 className="text-sm font-medium text-slate-400 mb-3 flex items-center gap-2">
                                        <FileText className="w-4 h-4" />
                                        Raw Data / Context
                                    </h4>
                                    <div className="text-sm text-slate-300 whitespace-pre-wrap font-mono leading-relaxed">
                                        {/* Mocking Raw Data for now since API might not return it fully yet */}
                                        {selectedDraft.raw_data || "Raw context data will appear here..."}
                                        {"\n\n"}
                                        [Excerpt from Source]
                                        {"\n"}
                                        The antibody-drug conjugate {selectedDraft.drug_name} targets {selectedDraft.target} and has shown promising results in early clinical trials...
                                    </div>
                                </div>

                                {/* Extracted JSON View */}
                                <div className="p-4 overflow-auto">
                                    <h4 className="text-sm font-medium text-purple-400 mb-3 flex items-center gap-2">
                                        <CheckCircle className="w-4 h-4" />
                                        Extracted Structured Data
                                    </h4>
                                    <pre className="text-sm text-slate-300 font-mono bg-slate-950 p-4 rounded-lg border border-slate-800 overflow-auto max-h-64">
                                        {JSON.stringify({
                                            name: selectedDraft.drug_name,
                                            target: selectedDraft.target,
                                            payload: selectedDraft.payload,
                                            linker: selectedDraft.linker,
                                            smiles_code: selectedDraft.smiles_code || null,
                                            outcome_type: selectedDraft.outcome_type,
                                            source: selectedDraft.enrichment_source,
                                            ai_analysis: selectedDraft.properties?.ai_analysis || null
                                        }, null, 2)}
                                    </pre>
                                </div>
                            </div>
                        </CardContent>
                    </>
                ) : (
                    <div className="flex-1 flex flex-col items-center justify-center text-slate-500">
                        <FileText className="w-12 h-12 mb-4 opacity-20" />
                        <p>Select an item from the inbox to review</p>
                    </div>
                )}
            </Card>

            {/* Reject Dialog */}
            <Dialog open={isRejectDialogOpen} onOpenChange={setIsRejectDialogOpen}>
                <DialogContent className="bg-slate-900 border-slate-800 text-white">
                    <DialogHeader>
                        <DialogTitle>Reject Draft</DialogTitle>
                        <DialogDescription className="text-slate-400">
                            Please provide a reason for rejecting this item. This will be logged.
                        </DialogDescription>
                    </DialogHeader>
                    <div className="py-4">
                        <Textarea
                            placeholder="Reason for rejection..."
                            className="bg-slate-800 border-slate-700 text-white min-h-[100px]"
                            value={rejectReason}
                            onChange={(e) => setRejectReason(e.target.value)}
                        />
                    </div>
                    <DialogFooter>
                        <Button
                            variant="outline"
                            onClick={() => setIsRejectDialogOpen(false)}
                            className="border-slate-700 text-slate-300"
                        >
                            Cancel
                        </Button>
                        <Button
                            variant="destructive"
                            onClick={() => selectedDraft && rejectMutation.mutate({ id: selectedDraft.id, reason: rejectReason })}
                            disabled={!rejectReason.trim()}
                        >
                            Confirm Reject
                        </Button>
                    </DialogFooter>
                </DialogContent>
            </Dialog>
        </div>
    )
}
