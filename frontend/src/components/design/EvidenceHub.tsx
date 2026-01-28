/**
 * Evidence Hub - Librarian 근거 논문 패널
 * 스캐폴드 클릭 시 관련 논문 문장 하이라이트
 */
import { useState, useEffect } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import { Button } from '@/components/ui/button'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Separator } from '@/components/ui/separator'
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger
} from '@/components/ui/tooltip'
import { motion, AnimatePresence } from 'framer-motion'
import {
  BookOpen,
  FileText,
  ExternalLink,
  ChevronRight,
  FlaskConical,
  Search,
  Sparkles,
  Link2
} from 'lucide-react'

export interface LiteratureReference {
  id: string
  pmid?: string
  title: string
  authors?: string
  journal?: string
  year?: number
  relevanceScore: number
  summary?: string
  highlightedSentences?: HighlightedSentence[]
  linkedScaffoldIds?: string[]
}

export interface HighlightedSentence {
  text: string
  relevance: 'high' | 'medium' | 'low'
  scaffoldId?: string
}

export interface GoldenSetReference {
  id: string
  drugName: string
  target: string
  indication?: string
  clinicalStatus: string
  similarity: number
}

interface EvidenceHubProps {
  references: LiteratureReference[]
  goldenSetRefs: GoldenSetReference[]
  selectedScaffoldId?: string | null
  isLoading?: boolean
  onReferenceClick?: (ref: LiteratureReference) => void
  onGoldenSetClick?: (ref: GoldenSetReference) => void
}

export function EvidenceHub({
  references,
  goldenSetRefs,
  selectedScaffoldId,
  isLoading = false,
  onReferenceClick,
  onGoldenSetClick
}: EvidenceHubProps) {
  const [activeTab, setActiveTab] = useState<'literature' | 'golden'>('literature')
  const [expandedRefId, setExpandedRefId] = useState<string | null>(null)

  // Filter references linked to selected scaffold
  const filteredRefs = selectedScaffoldId
    ? references.filter(r => r.linkedScaffoldIds?.includes(selectedScaffoldId))
    : references

  return (
    <Card className="h-full flex flex-col">
      <CardHeader className="pb-2">
        <div className="flex items-center justify-between">
          <CardTitle className="text-sm font-medium flex items-center gap-2">
            <BookOpen className="w-4 h-4 text-amber-500" />
            Evidence Hub
          </CardTitle>
          {selectedScaffoldId && (
            <Badge variant="outline" className="text-xs text-amber-600">
              <Link2 className="w-3 h-3 mr-1" />
              Linked to scaffold
            </Badge>
          )}
        </div>

        {/* Tab Navigation */}
        <div className="flex gap-1 mt-2">
          <Button
            variant={activeTab === 'literature' ? 'default' : 'ghost'}
            size="sm"
            className="h-7 text-xs"
            onClick={() => setActiveTab('literature')}
          >
            <FileText className="w-3 h-3 mr-1" />
            Literature ({filteredRefs.length})
          </Button>
          <Button
            variant={activeTab === 'golden' ? 'default' : 'ghost'}
            size="sm"
            className="h-7 text-xs"
            onClick={() => setActiveTab('golden')}
          >
            <FlaskConical className="w-3 h-3 mr-1" />
            Golden Set ({goldenSetRefs.length})
          </Button>
        </div>
      </CardHeader>

      <CardContent className="flex-1 overflow-hidden p-0">
        <ScrollArea className="h-[400px] p-4">
          {isLoading ? (
            <div className="flex flex-col items-center justify-center h-full text-gray-400">
              <Search className="w-8 h-8 animate-pulse mb-2" />
              <p className="text-sm">Searching knowledge base...</p>
            </div>
          ) : activeTab === 'literature' ? (
            <LiteratureList
              references={filteredRefs}
              expandedRefId={expandedRefId}
              selectedScaffoldId={selectedScaffoldId}
              onExpand={setExpandedRefId}
              onClick={onReferenceClick}
            />
          ) : (
            <GoldenSetList
              references={goldenSetRefs}
              onClick={onGoldenSetClick}
            />
          )}
        </ScrollArea>
      </CardContent>
    </Card>
  )
}

// Literature References List
function LiteratureList({
  references,
  expandedRefId,
  selectedScaffoldId,
  onExpand,
  onClick
}: {
  references: LiteratureReference[]
  expandedRefId: string | null
  selectedScaffoldId?: string | null
  onExpand: (id: string | null) => void
  onClick?: (ref: LiteratureReference) => void
}) {
  if (references.length === 0) {
    return (
      <div className="text-center text-gray-400 py-8">
        <FileText className="w-8 h-8 mx-auto mb-2 opacity-50" />
        <p className="text-sm">No literature references found</p>
      </div>
    )
  }

  return (
    <div className="space-y-2">
      <AnimatePresence>
        {references.map((ref, index) => (
          <motion.div
            key={ref.id}
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ delay: index * 0.05 }}
            className={`p-3 rounded-lg border cursor-pointer transition-colors ${
              expandedRefId === ref.id
                ? 'bg-amber-50 border-amber-200'
                : 'bg-white border-slate-200 hover:border-amber-200'
            }`}
            onClick={() => {
              onExpand(expandedRefId === ref.id ? null : ref.id)
              onClick?.(ref)
            }}
          >
            <div className="flex items-start justify-between gap-2">
              <div className="flex-1 min-w-0">
                <div className="flex items-center gap-2">
                  <Badge
                    variant="outline"
                    className={`text-[10px] shrink-0 ${
                      ref.relevanceScore >= 0.8
                        ? 'text-green-600 border-green-200'
                        : ref.relevanceScore >= 0.5
                        ? 'text-amber-600 border-amber-200'
                        : 'text-gray-600 border-gray-200'
                    }`}
                  >
                    {(ref.relevanceScore * 100).toFixed(0)}%
                  </Badge>
                  {ref.pmid && (
                    <Badge variant="secondary" className="text-[10px]">
                      PMID: {ref.pmid}
                    </Badge>
                  )}
                </div>
                <h4 className="text-xs font-medium mt-1 line-clamp-2">
                  {ref.title}
                </h4>
                {ref.authors && (
                  <p className="text-[10px] text-gray-500 mt-0.5 truncate">
                    {ref.authors}
                  </p>
                )}
              </div>
              <ChevronRight
                className={`w-4 h-4 text-gray-400 shrink-0 transition-transform ${
                  expandedRefId === ref.id ? 'rotate-90' : ''
                }`}
              />
            </div>

            {/* Expanded content with highlighted sentences */}
            {expandedRefId === ref.id && ref.highlightedSentences && (
              <motion.div
                initial={{ opacity: 0, height: 0 }}
                animate={{ opacity: 1, height: 'auto' }}
                exit={{ opacity: 0, height: 0 }}
                className="mt-3 pt-3 border-t border-amber-200"
              >
                <p className="text-[10px] text-amber-600 font-medium mb-2 flex items-center gap-1">
                  <Sparkles className="w-3 h-3" />
                  Key Evidence (Pinpoint)
                </p>
                <div className="space-y-2">
                  {ref.highlightedSentences.map((sentence, idx) => (
                    <div
                      key={idx}
                      className={`p-2 rounded text-[11px] ${
                        sentence.scaffoldId === selectedScaffoldId
                          ? 'bg-amber-100 border-l-2 border-amber-500'
                          : sentence.relevance === 'high'
                          ? 'bg-green-50 border-l-2 border-green-400'
                          : sentence.relevance === 'medium'
                          ? 'bg-yellow-50 border-l-2 border-yellow-400'
                          : 'bg-gray-50 border-l-2 border-gray-300'
                      }`}
                    >
                      "{sentence.text}"
                    </div>
                  ))}
                </div>

                {ref.pmid && (
                  <a
                    href={`https://pubmed.ncbi.nlm.nih.gov/${ref.pmid}`}
                    target="_blank"
                    rel="noopener noreferrer"
                    className="flex items-center gap-1 text-[10px] text-blue-600 hover:underline mt-2"
                    onClick={(e) => e.stopPropagation()}
                  >
                    <ExternalLink className="w-3 h-3" />
                    View on PubMed
                  </a>
                )}
              </motion.div>
            )}
          </motion.div>
        ))}
      </AnimatePresence>
    </div>
  )
}

// Golden Set References List
function GoldenSetList({
  references,
  onClick
}: {
  references: GoldenSetReference[]
  onClick?: (ref: GoldenSetReference) => void
}) {
  if (references.length === 0) {
    return (
      <div className="text-center text-gray-400 py-8">
        <FlaskConical className="w-8 h-8 mx-auto mb-2 opacity-50" />
        <p className="text-sm">No Golden Set matches found</p>
      </div>
    )
  }

  const getStatusColor = (status: string) => {
    if (status.toLowerCase().includes('approved')) return 'text-green-600 bg-green-50'
    if (status.toLowerCase().includes('phase 3')) return 'text-blue-600 bg-blue-50'
    if (status.toLowerCase().includes('phase 2')) return 'text-amber-600 bg-amber-50'
    return 'text-gray-600 bg-gray-50'
  }

  return (
    <div className="space-y-2">
      {references.map((ref, index) => (
        <motion.div
          key={ref.id}
          initial={{ opacity: 0, x: -10 }}
          animate={{ opacity: 1, x: 0 }}
          transition={{ delay: index * 0.05 }}
          className="p-3 rounded-lg border bg-white border-slate-200 hover:border-purple-200 cursor-pointer transition-colors"
          onClick={() => onClick?.(ref)}
        >
          <div className="flex items-start justify-between">
            <div>
              <h4 className="text-sm font-semibold text-purple-700">
                {ref.drugName}
              </h4>
              <p className="text-xs text-gray-600">
                Target: {ref.target}
              </p>
              {ref.indication && (
                <p className="text-[10px] text-gray-500">
                  {ref.indication}
                </p>
              )}
            </div>
            <div className="text-right">
              <Badge className={`text-[10px] ${getStatusColor(ref.clinicalStatus)}`}>
                {ref.clinicalStatus}
              </Badge>
              <p className="text-[10px] text-gray-500 mt-1">
                Similarity: {(ref.similarity * 100).toFixed(0)}%
              </p>
            </div>
          </div>
        </motion.div>
      ))}
    </div>
  )
}
