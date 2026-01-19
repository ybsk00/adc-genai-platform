import { useState } from 'react'
import { motion } from 'framer-motion'
import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query'
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
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
import { Textarea } from '@/components/ui/textarea'
import { ScrollArea } from '@/components/ui/scroll-area'
import { CheckCircle, XCircle, FileText, Loader2, ArrowRight } from 'lucide-react'
import { toast } from 'sonner'
import { supabase } from '@/lib/supabase'

interface GoldenSetDraft {
    id: string
    drug_name: string
    target: string
    payload?: string
    linker?: string
    enrichment_source: string
    created_at: string
    // Add other fields as needed for comparison
    raw_data?: string // Mock field for now if not in API response
}

export function StagingAreaTab() {
    const queryClient = useQueryClient()
    const [selectedDraft, setSelectedDraft] = useState<GoldenSetDraft | null>(null)
    const [rejectReason, setRejectReason] = useState('')
    const [isRejectDialogOpen, setIsRejectDialogOpen] = useState(false)

    // Fetch Drafts
    const { data: drafts, isLoading } = useQuery({
        queryKey: ['goldenSetDrafts'],
        queryFn: async () => {
            // Using direct supabase call or API endpoint
            // Since we have an API endpoint, let's use fetch or supabase client if configured for API
            // For now, let's use the supabase client directly if the API is just a wrapper, 
            // BUT the plan says "Use actual Backend API". 
            // Let's assume we have a fetch wrapper or just use fetch.
            // However, authentication might be tricky with raw fetch if not handled globally.
            // Let's use supabase-js to call the edge function or just standard fetch if we have a token.
            // Given the setup, I'll use the supabase client to get the session and then fetch.

            const { data: { session } } = await supabase.auth.getSession()
            if (!session) throw new Error('No session')

            const response = await fetch(`${import.meta.env.VITE_API_URL || 'http://localhost:8000'}/api/admin/goldenset/drafts`, {
                headers: {
                    Authorization: `Bearer ${session.access_token}`
                }
            })
            if (!response.ok) throw new Error('Failed to fetch drafts')
            return response.json() as Promise<GoldenSetDraft[]>
        }
    })

    // Approve Mutation
    const approveMutation = useMutation({
        mutationFn: async (id: string) => {
            const { data: { session } } = await supabase.auth.getSession()
            if (!session) throw new Error('No session')

            const response = await fetch(`${import.meta.env.VITE_API_URL || 'http://localhost:8000'}/api/admin/goldenset/${id}/approve`, {
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
            const { data: { session } } = await supabase.auth.getSession()
            if (!session) throw new Error('No session')

            const response = await fetch(`${import.meta.env.VITE_API_URL || 'http://localhost:8000'}/api/admin/goldenset/${id}/reject`, {
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
                <CardHeader>
                    <CardTitle className="text-white flex items-center justify-between">
                        <span>Inbox</span>
                        <Badge variant="secondary" className="bg-purple-500/20 text-purple-300">
                            {drafts?.length || 0} Pending
                        </Badge>
                    </CardTitle>
                    <CardDescription className="text-slate-400">
                        AI가 추출한 검토 대기 항목
                    </CardDescription>
                </CardHeader>
                <CardContent className="flex-1 overflow-hidden p-0">
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
                                        <Badge variant="outline" className="text-xs border-slate-700 text-slate-500">
                                            {draft.enrichment_source}
                                        </Badge>
                                    </div>
                                ))}
                            </div>
                        )}
                    </ScrollArea>
                </CardContent>
            </Card>

            {/* Comparison View */}
            <Card className="bg-slate-900 border-slate-800 lg:col-span-2 flex flex-col">
                {selectedDraft ? (
                    <>
                        <CardHeader className="border-b border-slate-800 pb-4">
                            <div className="flex justify-between items-center">
                                <div>
                                    <CardTitle className="text-white">{selectedDraft.drug_name}</CardTitle>
                                    <CardDescription className="text-slate-400 flex items-center gap-2 mt-1">
                                        Source: {selectedDraft.enrichment_source}
                                        <ArrowRight className="w-3 h-3" />
                                        Target: {selectedDraft.target}
                                    </CardDescription>
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
                                    <pre className="text-sm text-slate-300 font-mono bg-slate-950 p-4 rounded-lg border border-slate-800">
                                        {JSON.stringify({
                                            name: selectedDraft.drug_name,
                                            target: selectedDraft.target,
                                            payload: selectedDraft.payload,
                                            linker: selectedDraft.linker,
                                            source: selectedDraft.enrichment_source,
                                            status: 'draft'
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
