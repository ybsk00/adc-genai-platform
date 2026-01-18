import { useState, useCallback } from 'react'
import { useDropzone } from 'react-dropzone'
import { motion, AnimatePresence } from 'framer-motion'
import { Upload, FileText, X } from 'lucide-react'
import { Card, CardContent } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Progress } from '@/components/ui/progress'
import { toast } from 'sonner'
import { uploadDocument } from '@/lib/api'
import { supabase } from '@/lib/supabase'

interface UploadStatus {
    id: string
    file: File
    progress: number
    status: 'uploading' | 'processing' | 'success' | 'error'
    result?: any
    error?: string
}

export function UploadPage() {
    const [uploads, setUploads] = useState<UploadStatus[]>([])

    const onDrop = useCallback(async (acceptedFiles: File[]) => {
        const newUploads = acceptedFiles.map(file => ({
            id: Math.random().toString(36).substring(7),
            file,
            progress: 0,
            status: 'uploading' as const
        }))

        setUploads(prev => [...prev, ...newUploads])

        // 각 파일 업로드 처리
        for (const upload of newUploads) {
            await processUpload(upload)
        }
    }, [])

    const processUpload = async (upload: UploadStatus) => {
        try {
            // 1. Uploading
            updateUploadStatus(upload.id, { status: 'uploading', progress: 30 })

            // Get Session Token
            const { data: { session } } = await supabase?.auth.getSession() || { data: { session: null } }
            const token = session?.access_token || ''

            if (!token) {
                throw new Error('로그인이 필요합니다.')
            }

            // 2. Processing (Gemini Parsing)
            updateUploadStatus(upload.id, { status: 'processing', progress: 60 })

            const { data, error } = await uploadDocument(upload.file, {}, token)

            if (error) {
                // Check for specific error messages if needed, e.g. duplicate
                if (error.includes('duplicate') || error.includes('중복')) {
                    updateUploadStatus(upload.id, {
                        status: 'error',
                        progress: 100,
                        error: `중복된 문서입니다: ${error}`
                    })
                    toast.warning(`${upload.file.name}: 중복된 문서`)
                    return
                }
                throw new Error(error)
            }

            if (data) {
                updateUploadStatus(upload.id, {
                    status: 'success',
                    progress: 100,
                    result: data
                })
                toast.success(`${upload.file.name} 파싱 완료`)
            }

        } catch (error: any) {
            updateUploadStatus(upload.id, {
                status: 'error',
                progress: 100,
                error: error.message || '업로드 실패'
            })
            toast.error(`${upload.file.name} 처리 실패`)
        }
    }

    const updateUploadStatus = (id: string, updates: Partial<UploadStatus>) => {
        setUploads(prev => prev.map(u => u.id === id ? { ...u, ...updates } : u))
    }

    const removeUpload = (id: string) => {
        setUploads(prev => prev.filter(u => u.id !== id))
    }

    const { getRootProps, getInputProps, isDragActive } = useDropzone({
        onDrop,
        accept: {
            'application/pdf': ['.pdf']
        },
        multiple: true
    })

    return (
        <div className="space-y-6">
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
            >
                <h1 className="text-2xl font-bold text-white">Document Upload</h1>
                <p className="text-slate-400 mt-1">
                    PDF 문서를 업로드하여 Gemini Vision으로 구조화된 데이터를 추출합니다.
                </p>
            </motion.div>

            {/* Dropzone */}
            <Card className="bg-slate-900 border-slate-800">
                <CardContent className="p-0">
                    <div
                        {...getRootProps()}
                        className={`
                            border-2 border-dashed rounded-lg p-12 text-center cursor-pointer transition-colors
                            ${isDragActive
                                ? 'border-purple-500 bg-purple-500/10'
                                : 'border-slate-700 hover:border-slate-600 hover:bg-slate-800/50'
                            }
                        `}
                    >
                        <input {...getInputProps()} />
                        <div className="flex flex-col items-center gap-4">
                            <div className="p-4 bg-slate-800 rounded-full">
                                <Upload className="w-8 h-8 text-slate-400" />
                            </div>
                            <div>
                                <p className="text-lg font-medium text-white">
                                    {isDragActive ? '파일을 놓으세요' : 'PDF 파일을 드래그하거나 클릭하여 업로드'}
                                </p>
                                <p className="text-sm text-slate-500 mt-1">
                                    최대 10MB, PDF 형식만 지원
                                </p>
                            </div>
                        </div>
                    </div>
                </CardContent>
            </Card>

            {/* Upload List */}
            <AnimatePresence>
                {uploads.length > 0 && (
                    <motion.div
                        initial={{ opacity: 0, height: 0 }}
                        animate={{ opacity: 1, height: 'auto' }}
                        exit={{ opacity: 0, height: 0 }}
                        className="space-y-3"
                    >
                        {uploads.map(upload => (
                            <motion.div
                                key={upload.id}
                                initial={{ opacity: 0, x: -20 }}
                                animate={{ opacity: 1, x: 0 }}
                                exit={{ opacity: 0, x: 20 }}
                            >
                                <Card className="bg-slate-900 border-slate-800">
                                    <CardContent className="p-4 flex items-center gap-4">
                                        <div className="p-2 bg-slate-800 rounded-lg">
                                            <FileText className="w-6 h-6 text-slate-400" />
                                        </div>

                                        <div className="flex-1 min-w-0">
                                            <div className="flex justify-between mb-1">
                                                <p className="font-medium text-white truncate">
                                                    {upload.file.name}
                                                </p>
                                                <span className={`text-xs ${upload.status === 'error' ? 'text-red-400' :
                                                    upload.status === 'success' ? 'text-green-400' :
                                                        'text-blue-400'
                                                    }`}>
                                                    {upload.status === 'uploading' && 'Uploading...'}
                                                    {upload.status === 'processing' && 'Gemini Parsing...'}
                                                    {upload.status === 'success' && 'Completed'}
                                                    {upload.status === 'error' && 'Failed'}
                                                </span>
                                            </div>
                                            <Progress
                                                value={upload.progress}
                                                className={`h-1.5 ${upload.status === 'error' ? 'bg-red-900' : 'bg-slate-800'
                                                    }`}
                                            />
                                            {upload.error && (
                                                <p className="text-xs text-red-400 mt-1">
                                                    {upload.error}
                                                </p>
                                            )}
                                        </div>

                                        <Button
                                            variant="ghost"
                                            size="icon"
                                            className="text-slate-500 hover:text-white"
                                            onClick={() => removeUpload(upload.id)}
                                        >
                                            <X className="w-4 h-4" />
                                        </Button>
                                    </CardContent>
                                </Card>
                            </motion.div>
                        ))}
                    </motion.div>
                )}
            </AnimatePresence>
        </div>
    )
}
