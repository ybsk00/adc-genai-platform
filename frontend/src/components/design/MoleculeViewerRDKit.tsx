/**
 * MoleculeViewerRDKit Component
 * RDKit.js WASM 기반 고품질 분자 구조 뷰어
 * smiles-drawer보다 정확한 화학 구조 렌더링 제공
 */
import { useState, useEffect, useRef, useCallback } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import { Button } from '@/components/ui/button'
import { Skeleton } from '@/components/ui/skeleton'
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger
} from '@/components/ui/tooltip'
import {
  FlaskConical,
  Download,
  Copy,
  ZoomIn,
  ZoomOut,
  RotateCcw,
  Loader2,
  AlertCircle
} from 'lucide-react'
import { toast } from 'sonner'
import { cn } from '@/lib/utils'

// RDKit.js 타입 정의
interface RDKitModule {
  get_mol: (smiles: string) => RDKitMolecule | null
  prefer_coordgen: (prefer: boolean) => void
  version: () => string
}

interface RDKitMolecule {
  is_valid: () => boolean
  get_svg: (width?: number, height?: number) => string
  get_svg_with_highlights: (details: string) => string
  get_molblock: () => string
  get_inchi: () => string
  get_descriptors: () => string
  delete: () => void
}

// 전역 RDKit 인스턴스
let rdkitInstance: RDKitModule | null = null
let rdkitLoadPromise: Promise<RDKitModule> | null = null

// RDKit.js 로드 함수
async function loadRDKit(): Promise<RDKitModule> {
  if (rdkitInstance) return rdkitInstance

  if (rdkitLoadPromise) return rdkitLoadPromise

  rdkitLoadPromise = new Promise((resolve, reject) => {
    // RDKit.js CDN에서 로드
    const script = document.createElement('script')
    script.src = 'https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js'
    script.async = true

    script.onload = async () => {
      try {
        // @ts-ignore - RDKit global
        const initRDKit = window.initRDKitModule
        if (!initRDKit) {
          throw new Error('RDKit module not found')
        }

        const rdkit = await initRDKit()
        rdkit.prefer_coordgen(true) // 더 나은 좌표 생성
        rdkitInstance = rdkit
        console.log('[RDKit] Loaded version:', rdkit.version())
        resolve(rdkit)
      } catch (err) {
        reject(err)
      }
    }

    script.onerror = () => reject(new Error('Failed to load RDKit.js'))
    document.head.appendChild(script)
  })

  return rdkitLoadPromise
}

export interface MoleculeViewerRDKitProps {
  smiles: string
  title?: string
  width?: number
  height?: number
  showControls?: boolean
  showDescriptors?: boolean
  highlightAtoms?: number[]
  highlightBonds?: number[]
  highlightColor?: string
  className?: string
  onLoad?: (descriptors: MoleculeDescriptors) => void
  onError?: (error: string) => void
}

export interface MoleculeDescriptors {
  mw: number
  logP: number
  hbd: number
  hba: number
  tpsa: number
  rotatable_bonds: number
  num_rings: number
  num_aromatic_rings: number
  fraction_csp3: number
}

export function MoleculeViewerRDKit({
  smiles,
  title,
  width = 400,
  height = 300,
  showControls = true,
  showDescriptors = false,
  highlightAtoms = [],
  highlightBonds = [],
  highlightColor = '#3b82f6',
  className,
  onLoad,
  onError
}: MoleculeViewerRDKitProps) {
  const [svgContent, setSvgContent] = useState<string>('')
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const [descriptors, setDescriptors] = useState<MoleculeDescriptors | null>(null)
  const [zoom, setZoom] = useState(1)
  const containerRef = useRef<HTMLDivElement>(null)

  // 분자 렌더링
  const renderMolecule = useCallback(async () => {
    if (!smiles) {
      setSvgContent('')
      setIsLoading(false)
      return
    }

    setIsLoading(true)
    setError(null)

    try {
      const rdkit = await loadRDKit()
      const mol = rdkit.get_mol(smiles)

      if (!mol || !mol.is_valid()) {
        throw new Error('Invalid SMILES structure')
      }

      // 하이라이트 옵션 생성
      let svg: string
      if (highlightAtoms.length > 0 || highlightBonds.length > 0) {
        const highlightDetails = JSON.stringify({
          atoms: highlightAtoms,
          bonds: highlightBonds,
          highlightColour: hexToRgb(highlightColor),
          width: width * zoom,
          height: height * zoom
        })
        svg = mol.get_svg_with_highlights(highlightDetails)
      } else {
        svg = mol.get_svg(width * zoom, height * zoom)
      }

      setSvgContent(svg)

      // 분자 설명자 계산
      try {
        const descJson = mol.get_descriptors()
        const desc = JSON.parse(descJson)
        const moleculeDescriptors: MoleculeDescriptors = {
          mw: desc.exactmw || desc.amw || 0,
          logP: desc.CrippenClogP || 0,
          hbd: desc.NumHBD || 0,
          hba: desc.NumHBA || 0,
          tpsa: desc.TPSA || 0,
          rotatable_bonds: desc.NumRotatableBonds || 0,
          num_rings: desc.NumRings || 0,
          num_aromatic_rings: desc.NumAromaticRings || 0,
          fraction_csp3: desc.FractionCSP3 || 0
        }
        setDescriptors(moleculeDescriptors)
        onLoad?.(moleculeDescriptors)
      } catch {
        // 설명자 계산 실패 시 무시
      }

      mol.delete() // 메모리 해제
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'Failed to render molecule'
      setError(errorMessage)
      onError?.(errorMessage)
    } finally {
      setIsLoading(false)
    }
  }, [smiles, width, height, zoom, highlightAtoms, highlightBonds, highlightColor, onLoad, onError])

  useEffect(() => {
    renderMolecule()
  }, [renderMolecule])

  // 줌 핸들러
  const handleZoomIn = () => setZoom(prev => Math.min(prev + 0.25, 3))
  const handleZoomOut = () => setZoom(prev => Math.max(prev - 0.25, 0.5))
  const handleZoomReset = () => setZoom(1)

  // 복사 핸들러
  const handleCopySMILES = async () => {
    await navigator.clipboard.writeText(smiles)
    toast.success('SMILES copied to clipboard')
  }

  // SVG 다운로드 핸들러
  const handleDownloadSVG = () => {
    if (!svgContent) return

    const blob = new Blob([svgContent], { type: 'image/svg+xml' })
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = `molecule-${Date.now()}.svg`
    a.click()
    URL.revokeObjectURL(url)
    toast.success('SVG downloaded')
  }

  return (
    <Card className={cn('overflow-hidden', className)}>
      {title && (
        <CardHeader className="pb-2">
          <CardTitle className="text-sm flex items-center gap-2">
            <FlaskConical className="w-4 h-4" />
            {title}
            <Badge variant="outline" className="ml-auto text-xs">
              RDKit.js
            </Badge>
          </CardTitle>
        </CardHeader>
      )}

      <CardContent className="p-0">
        {/* 컨트롤 바 */}
        {showControls && (
          <div className="flex items-center justify-between px-3 py-2 border-b bg-gray-50">
            <div className="flex items-center gap-1">
              <TooltipProvider>
                <Tooltip>
                  <TooltipTrigger asChild>
                    <Button variant="ghost" size="icon" className="h-7 w-7" onClick={handleZoomIn}>
                      <ZoomIn className="w-4 h-4" />
                    </Button>
                  </TooltipTrigger>
                  <TooltipContent>Zoom In</TooltipContent>
                </Tooltip>

                <Tooltip>
                  <TooltipTrigger asChild>
                    <Button variant="ghost" size="icon" className="h-7 w-7" onClick={handleZoomOut}>
                      <ZoomOut className="w-4 h-4" />
                    </Button>
                  </TooltipTrigger>
                  <TooltipContent>Zoom Out</TooltipContent>
                </Tooltip>

                <Tooltip>
                  <TooltipTrigger asChild>
                    <Button variant="ghost" size="icon" className="h-7 w-7" onClick={handleZoomReset}>
                      <RotateCcw className="w-4 h-4" />
                    </Button>
                  </TooltipTrigger>
                  <TooltipContent>Reset</TooltipContent>
                </Tooltip>
              </TooltipProvider>

              <span className="text-xs text-muted-foreground ml-2">
                {Math.round(zoom * 100)}%
              </span>
            </div>

            <div className="flex items-center gap-1">
              <TooltipProvider>
                <Tooltip>
                  <TooltipTrigger asChild>
                    <Button variant="ghost" size="icon" className="h-7 w-7" onClick={handleCopySMILES}>
                      <Copy className="w-4 h-4" />
                    </Button>
                  </TooltipTrigger>
                  <TooltipContent>Copy SMILES</TooltipContent>
                </Tooltip>

                <Tooltip>
                  <TooltipTrigger asChild>
                    <Button variant="ghost" size="icon" className="h-7 w-7" onClick={handleDownloadSVG}>
                      <Download className="w-4 h-4" />
                    </Button>
                  </TooltipTrigger>
                  <TooltipContent>Download SVG</TooltipContent>
                </Tooltip>
              </TooltipProvider>
            </div>
          </div>
        )}

        {/* 분자 뷰어 */}
        <div
          ref={containerRef}
          className="flex items-center justify-center bg-white overflow-auto"
          style={{ minHeight: height }}
        >
          {isLoading ? (
            <div className="flex flex-col items-center gap-2">
              <Loader2 className="w-8 h-8 animate-spin text-blue-500" />
              <span className="text-sm text-muted-foreground">Loading RDKit.js...</span>
            </div>
          ) : error ? (
            <div className="flex flex-col items-center gap-2 text-red-500">
              <AlertCircle className="w-8 h-8" />
              <span className="text-sm">{error}</span>
            </div>
          ) : svgContent ? (
            <div
              dangerouslySetInnerHTML={{ __html: svgContent }}
              className="molecule-svg"
            />
          ) : (
            <div className="text-center text-muted-foreground">
              <FlaskConical className="w-12 h-12 mx-auto mb-2 opacity-30" />
              <span className="text-sm">No structure to display</span>
            </div>
          )}
        </div>

        {/* 설명자 표시 */}
        {showDescriptors && descriptors && (
          <div className="grid grid-cols-4 gap-2 p-3 border-t bg-gray-50 text-center">
            <DescriptorItem label="MW" value={descriptors.mw.toFixed(1)} />
            <DescriptorItem label="LogP" value={descriptors.logP.toFixed(2)} />
            <DescriptorItem label="HBD" value={descriptors.hbd.toString()} />
            <DescriptorItem label="HBA" value={descriptors.hba.toString()} />
            <DescriptorItem label="TPSA" value={descriptors.tpsa.toFixed(1)} />
            <DescriptorItem label="RotB" value={descriptors.rotatable_bonds.toString()} />
            <DescriptorItem label="Rings" value={descriptors.num_rings.toString()} />
            <DescriptorItem label="Fsp3" value={descriptors.fraction_csp3.toFixed(2)} />
          </div>
        )}
      </CardContent>
    </Card>
  )
}

// 설명자 아이템 컴포넌트
function DescriptorItem({ label, value }: { label: string; value: string }) {
  return (
    <div className="text-xs">
      <span className="text-muted-foreground">{label}</span>
      <div className="font-mono font-medium">{value}</div>
    </div>
  )
}

// 헥스 색상을 RGB 배열로 변환
function hexToRgb(hex: string): [number, number, number] {
  const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex)
  return result
    ? [parseInt(result[1], 16) / 255, parseInt(result[2], 16) / 255, parseInt(result[3], 16) / 255]
    : [0.23, 0.51, 0.96] // 기본 파란색
}

// RDKit 사용 가능 여부 확인 훅
export function useRDKitAvailable(): boolean {
  const [available, setAvailable] = useState(false)

  useEffect(() => {
    loadRDKit()
      .then(() => setAvailable(true))
      .catch(() => setAvailable(false))
  }, [])

  return available
}

export default MoleculeViewerRDKit
