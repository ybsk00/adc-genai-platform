/**
 * AntigenAutocomplete Component
 * UniProt API를 활용한 Target Antigen 자동완성 컴포넌트
 */
import { useState, useEffect, useRef, useCallback } from 'react'
import { Input } from '@/components/ui/input'
import { Badge } from '@/components/ui/badge'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Loader2, Search, ExternalLink, Dna } from 'lucide-react'
import { motion, AnimatePresence } from 'framer-motion'
import { cn } from '@/lib/utils'

const API_BASE_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000'

interface AntigenResult {
  uniprot_id: string
  name: string
  gene_name: string | null
  organism: string
  sequence_length?: number
}

interface AntigenAutocompleteProps {
  value: string
  onChange: (value: string, selected?: AntigenResult) => void
  placeholder?: string
  disabled?: boolean
  className?: string
}

// 디바운스 훅
function useDebounce<T>(value: T, delay: number): T {
  const [debouncedValue, setDebouncedValue] = useState<T>(value)

  useEffect(() => {
    const timer = setTimeout(() => setDebouncedValue(value), delay)
    return () => clearTimeout(timer)
  }, [value, delay])

  return debouncedValue
}

export function AntigenAutocomplete({
  value,
  onChange,
  placeholder = 'Search antigen (e.g., HER2, TROP2, EGFR)',
  disabled = false,
  className
}: AntigenAutocompleteProps) {
  const [query, setQuery] = useState(value)
  const [results, setResults] = useState<AntigenResult[]>([])
  const [commonTargets, setCommonTargets] = useState<AntigenResult[]>([])
  const [isLoading, setIsLoading] = useState(false)
  const [isOpen, setIsOpen] = useState(false)
  const [selectedIndex, setSelectedIndex] = useState(-1)
  const inputRef = useRef<HTMLInputElement>(null)
  const containerRef = useRef<HTMLDivElement>(null)

  const debouncedQuery = useDebounce(query, 300)

  // 자주 사용되는 타겟 로드
  useEffect(() => {
    const loadCommonTargets = async () => {
      try {
        const response = await fetch(`${API_BASE_URL}/api/uniprot/common-targets`)
        if (response.ok) {
          const data = await response.json()
          setCommonTargets(data.targets || [])
        }
      } catch (error) {
        console.error('Failed to load common targets:', error)
      }
    }
    loadCommonTargets()
  }, [])

  // UniProt 검색
  useEffect(() => {
    const searchAntigens = async () => {
      if (debouncedQuery.length < 2) {
        setResults([])
        return
      }

      setIsLoading(true)
      try {
        const response = await fetch(
          `${API_BASE_URL}/api/uniprot/search?query=${encodeURIComponent(debouncedQuery)}&limit=10`
        )
        if (response.ok) {
          const data = await response.json()
          setResults(data.results || [])
        }
      } catch (error) {
        console.error('UniProt search failed:', error)
        // Fallback to common targets filter
        const filtered = commonTargets.filter(t =>
          t.gene_name?.toLowerCase().includes(debouncedQuery.toLowerCase()) ||
          t.name.toLowerCase().includes(debouncedQuery.toLowerCase())
        )
        setResults(filtered)
      } finally {
        setIsLoading(false)
      }
    }

    searchAntigens()
  }, [debouncedQuery, commonTargets])

  // 외부 클릭 감지
  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (containerRef.current && !containerRef.current.contains(event.target as Node)) {
        setIsOpen(false)
      }
    }

    document.addEventListener('mousedown', handleClickOutside)
    return () => document.removeEventListener('mousedown', handleClickOutside)
  }, [])

  // 키보드 네비게이션
  const handleKeyDown = useCallback((e: React.KeyboardEvent) => {
    const items = results.length > 0 ? results : commonTargets.slice(0, 5)

    switch (e.key) {
      case 'ArrowDown':
        e.preventDefault()
        setSelectedIndex(prev => Math.min(prev + 1, items.length - 1))
        break
      case 'ArrowUp':
        e.preventDefault()
        setSelectedIndex(prev => Math.max(prev - 1, -1))
        break
      case 'Enter':
        e.preventDefault()
        if (selectedIndex >= 0 && items[selectedIndex]) {
          handleSelect(items[selectedIndex])
        }
        break
      case 'Escape':
        setIsOpen(false)
        break
    }
  }, [results, commonTargets, selectedIndex])

  const handleSelect = (item: AntigenResult) => {
    const displayValue = item.gene_name || item.name
    setQuery(displayValue)
    onChange(displayValue, item)
    setIsOpen(false)
    setSelectedIndex(-1)
  }

  const handleInputChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = e.target.value
    setQuery(newValue)
    onChange(newValue)
    setIsOpen(true)
    setSelectedIndex(-1)
  }

  const displayItems = results.length > 0 ? results : (query.length === 0 ? commonTargets.slice(0, 8) : [])

  return (
    <div ref={containerRef} className={cn('relative', className)}>
      <div className="relative">
        <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-muted-foreground" />
        <Input
          ref={inputRef}
          value={query}
          onChange={handleInputChange}
          onFocus={() => setIsOpen(true)}
          onKeyDown={handleKeyDown}
          placeholder={placeholder}
          disabled={disabled}
          className="pl-10 pr-10"
        />
        {isLoading && (
          <Loader2 className="absolute right-3 top-1/2 -translate-y-1/2 w-4 h-4 animate-spin text-muted-foreground" />
        )}
      </div>

      <AnimatePresence>
        {isOpen && displayItems.length > 0 && (
          <motion.div
            initial={{ opacity: 0, y: -10 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -10 }}
            transition={{ duration: 0.15 }}
            className="absolute z-50 w-full mt-1 bg-white border rounded-lg shadow-lg overflow-hidden"
          >
            {query.length === 0 && (
              <div className="px-3 py-2 text-xs text-muted-foreground bg-gray-50 border-b">
                Common ADC Targets
              </div>
            )}
            <ScrollArea className="max-h-64">
              {displayItems.map((item, index) => (
                <div
                  key={`${item.uniprot_id}-${index}`}
                  onClick={() => handleSelect(item)}
                  className={cn(
                    'px-3 py-2 cursor-pointer transition-colors',
                    selectedIndex === index
                      ? 'bg-blue-50 border-l-2 border-blue-500'
                      : 'hover:bg-gray-50'
                  )}
                >
                  <div className="flex items-center justify-between">
                    <div className="flex items-center gap-2">
                      <Dna className="w-4 h-4 text-blue-500" />
                      <span className="font-medium">
                        {item.gene_name || item.name}
                      </span>
                      {item.gene_name && (
                        <span className="text-xs text-muted-foreground truncate max-w-[150px]">
                          {item.name}
                        </span>
                      )}
                    </div>
                    <Badge variant="outline" className="text-xs shrink-0">
                      {item.uniprot_id}
                    </Badge>
                  </div>
                  <div className="flex items-center gap-2 mt-1 text-xs text-muted-foreground">
                    <span>{item.organism}</span>
                    {item.sequence_length && (
                      <span>| {item.sequence_length} aa</span>
                    )}
                  </div>
                </div>
              ))}
            </ScrollArea>
            <div className="px-3 py-2 text-xs text-muted-foreground bg-gray-50 border-t flex items-center gap-1">
              <span>Data from</span>
              <a
                href="https://www.uniprot.org"
                target="_blank"
                rel="noopener noreferrer"
                className="text-blue-500 hover:underline flex items-center gap-1"
              >
                UniProt
                <ExternalLink className="w-3 h-3" />
              </a>
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  )
}

// Compact 버전 (인라인 표시용)
export function AntigenBadge({ antigen, uniprotId }: { antigen: string; uniprotId?: string }) {
  return (
    <Badge variant="secondary" className="gap-1">
      <Dna className="w-3 h-3" />
      {antigen}
      {uniprotId && (
        <a
          href={`https://www.uniprot.org/uniprotkb/${uniprotId}`}
          target="_blank"
          rel="noopener noreferrer"
          className="ml-1 text-blue-500 hover:underline"
          onClick={(e) => e.stopPropagation()}
        >
          <ExternalLink className="w-3 h-3" />
        </a>
      )}
    </Badge>
  )
}

export default AntigenAutocomplete
