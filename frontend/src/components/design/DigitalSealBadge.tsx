/**
 * Digital Seal Badge - FDA 21 CFR Part 11 준수 해시 배지
 * 데이터 무결성 검증 시각화
 */
import { useState } from 'react'
import { Card, CardContent } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import { Button } from '@/components/ui/button'
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger
} from '@/components/ui/tooltip'
import {
  Dialog,
  DialogContent,
  DialogDescription,
  DialogHeader,
  DialogTitle,
  DialogTrigger
} from '@/components/ui/dialog'
import { motion } from 'framer-motion'
import {
  Shield,
  ShieldCheck,
  ShieldAlert,
  Copy,
  Check,
  FileText,
  Download,
  ExternalLink,
  Clock,
  Hash,
  Link2,
  Lock
} from 'lucide-react'

interface DigitalSealBadgeProps {
  sessionId: string
  recordHash: string
  previousHash?: string
  chainHash?: string
  timestamp: string
  isVerified?: boolean
  onExportBundle?: () => void
  className?: string
}

export function DigitalSealBadge({
  sessionId,
  recordHash,
  previousHash,
  chainHash,
  timestamp,
  isVerified = true,
  onExportBundle,
  className
}: DigitalSealBadgeProps) {
  const [copied, setCopied] = useState(false)
  const [showDetails, setShowDetails] = useState(false)

  const handleCopy = async (text: string) => {
    await navigator.clipboard.writeText(text)
    setCopied(true)
    setTimeout(() => setCopied(false), 2000)
  }

  const formatHash = (hash: string) => {
    if (!hash) return 'N/A'
    return `${hash.slice(0, 8)}...${hash.slice(-8)}`
  }

  const formatTimestamp = (ts: string) => {
    const date = new Date(ts)
    return date.toLocaleString('ko-KR', {
      year: 'numeric',
      month: '2-digit',
      day: '2-digit',
      hour: '2-digit',
      minute: '2-digit',
      second: '2-digit'
    })
  }

  return (
    <Card className={className}>
      <CardContent className="p-3">
        <div className="flex items-center justify-between">
          {/* Badge */}
          <TooltipProvider>
            <Tooltip>
              <TooltipTrigger asChild>
                <motion.div
                  initial={{ scale: 0.9, opacity: 0 }}
                  animate={{ scale: 1, opacity: 1 }}
                  className={`flex items-center gap-2 px-3 py-1.5 rounded-full border ${
                    isVerified
                      ? 'bg-green-50 border-green-200 text-green-700'
                      : 'bg-red-50 border-red-200 text-red-700'
                  }`}
                >
                  {isVerified ? (
                    <ShieldCheck className="w-4 h-4" />
                  ) : (
                    <ShieldAlert className="w-4 h-4" />
                  )}
                  <span className="text-xs font-medium">
                    21 CFR Part 11 Compliant
                  </span>
                </motion.div>
              </TooltipTrigger>
              <TooltipContent>
                <p>Data integrity verified via chain hash</p>
              </TooltipContent>
            </Tooltip>
          </TooltipProvider>

          {/* Hash Preview */}
          <div className="flex items-center gap-2">
            <Badge
              variant="outline"
              className="font-mono text-[10px] cursor-pointer hover:bg-slate-50"
              onClick={() => handleCopy(recordHash)}
            >
              <Hash className="w-3 h-3 mr-1" />
              {formatHash(recordHash)}
              {copied ? (
                <Check className="w-3 h-3 ml-1 text-green-500" />
              ) : (
                <Copy className="w-3 h-3 ml-1 opacity-50" />
              )}
            </Badge>

            {/* Details Dialog */}
            <Dialog open={showDetails} onOpenChange={setShowDetails}>
              <DialogTrigger asChild>
                <Button variant="ghost" size="icon" className="h-7 w-7">
                  <FileText className="w-3 h-3" />
                </Button>
              </DialogTrigger>
              <DialogContent className="max-w-md">
                <DialogHeader>
                  <DialogTitle className="flex items-center gap-2">
                    <Shield className="w-5 h-5 text-green-600" />
                    Digital Seal Details
                  </DialogTitle>
                  <DialogDescription>
                    FDA 21 CFR Part 11 Compliance Record
                  </DialogDescription>
                </DialogHeader>

                <div className="space-y-4 mt-4">
                  {/* Session ID */}
                  <div className="p-3 bg-slate-50 rounded-lg">
                    <p className="text-[10px] text-gray-500 uppercase mb-1">
                      Session ID
                    </p>
                    <p className="font-mono text-xs break-all">{sessionId}</p>
                  </div>

                  {/* Record Hash */}
                  <div className="p-3 bg-green-50 rounded-lg border border-green-200">
                    <p className="text-[10px] text-green-600 uppercase mb-1 flex items-center gap-1">
                      <Lock className="w-3 h-3" />
                      Record Hash (SHA-256)
                    </p>
                    <p className="font-mono text-xs break-all text-green-800">
                      {recordHash}
                    </p>
                  </div>

                  {/* Previous Hash */}
                  {previousHash && (
                    <div className="p-3 bg-slate-50 rounded-lg">
                      <p className="text-[10px] text-gray-500 uppercase mb-1 flex items-center gap-1">
                        <Link2 className="w-3 h-3" />
                        Previous Hash (Chain Link)
                      </p>
                      <p className="font-mono text-xs break-all">
                        {previousHash}
                      </p>
                    </div>
                  )}

                  {/* Chain Hash */}
                  {chainHash && (
                    <div className="p-3 bg-blue-50 rounded-lg border border-blue-200">
                      <p className="text-[10px] text-blue-600 uppercase mb-1">
                        Chain Hash
                      </p>
                      <p className="font-mono text-xs break-all text-blue-800">
                        {chainHash}
                      </p>
                    </div>
                  )}

                  {/* Timestamp */}
                  <div className="p-3 bg-slate-50 rounded-lg">
                    <p className="text-[10px] text-gray-500 uppercase mb-1 flex items-center gap-1">
                      <Clock className="w-3 h-3" />
                      Timestamp (UTC)
                    </p>
                    <p className="font-mono text-xs">{formatTimestamp(timestamp)}</p>
                  </div>

                  {/* Verification Status */}
                  <div className={`p-3 rounded-lg ${
                    isVerified
                      ? 'bg-green-100 border border-green-300'
                      : 'bg-red-100 border border-red-300'
                  }`}>
                    <div className="flex items-center gap-2">
                      {isVerified ? (
                        <>
                          <ShieldCheck className="w-5 h-5 text-green-600" />
                          <div>
                            <p className="text-sm font-medium text-green-800">
                              Integrity Verified
                            </p>
                            <p className="text-[10px] text-green-600">
                              Hash chain validated successfully
                            </p>
                          </div>
                        </>
                      ) : (
                        <>
                          <ShieldAlert className="w-5 h-5 text-red-600" />
                          <div>
                            <p className="text-sm font-medium text-red-800">
                              Verification Failed
                            </p>
                            <p className="text-[10px] text-red-600">
                              Hash mismatch detected
                            </p>
                          </div>
                        </>
                      )}
                    </div>
                  </div>

                  {/* Export Button */}
                  {onExportBundle && (
                    <Button
                      className="w-full"
                      variant="outline"
                      onClick={onExportBundle}
                    >
                      <Download className="w-4 h-4 mr-2" />
                      Download Audit Trail Bundle
                    </Button>
                  )}
                </div>
              </DialogContent>
            </Dialog>

            {/* Quick Export */}
            {onExportBundle && (
              <TooltipProvider>
                <Tooltip>
                  <TooltipTrigger asChild>
                    <Button
                      variant="ghost"
                      size="icon"
                      className="h-7 w-7"
                      onClick={onExportBundle}
                    >
                      <Download className="w-3 h-3" />
                    </Button>
                  </TooltipTrigger>
                  <TooltipContent>
                    <p>Download Audit Trail Bundle</p>
                  </TooltipContent>
                </Tooltip>
              </TooltipProvider>
            )}
          </div>
        </div>
      </CardContent>
    </Card>
  )
}

/**
 * Compact version for inline use
 */
interface CompactSealBadgeProps {
  hash: string
  isVerified?: boolean
  onClick?: () => void
}

export function CompactSealBadge({
  hash,
  isVerified = true,
  onClick
}: CompactSealBadgeProps) {
  return (
    <Badge
      variant="outline"
      className={`cursor-pointer ${
        isVerified
          ? 'text-green-600 border-green-200 hover:bg-green-50'
          : 'text-red-600 border-red-200 hover:bg-red-50'
      }`}
      onClick={onClick}
    >
      {isVerified ? (
        <ShieldCheck className="w-3 h-3 mr-1" />
      ) : (
        <ShieldAlert className="w-3 h-3 mr-1" />
      )}
      <span className="font-mono text-[10px]">
        {hash.slice(0, 8)}...
      </span>
    </Badge>
  )
}
