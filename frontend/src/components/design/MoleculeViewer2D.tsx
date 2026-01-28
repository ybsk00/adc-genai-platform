/**
 * Molecule Viewer 2D - RDKit.js WASM 기반 분자 구조 시각화
 * 업그레이드: smiles-drawer → RDKit.js WASM
 */
import { useEffect, useRef, useState, useCallback } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import { Skeleton } from '@/components/ui/skeleton'
import {
  Atom,
  ZoomIn,
  ZoomOut,
  RotateCcw,
  Download,
  Copy,
  Check,
  AlertCircle
} from 'lucide-react'
import type { RDKitModule, JSMol } from '@rdkit/rdkit'

// Global RDKit instance (singleton)
let rdkitInstance: RDKitModule | null = null
let rdkitLoadPromise: Promise<RDKitModule> | null = null

// Initialize RDKit WASM
async function initRDKit(): Promise<RDKitModule> {
  if (rdkitInstance) return rdkitInstance

  if (rdkitLoadPromise) return rdkitLoadPromise

  rdkitLoadPromise = (async () => {
    try {
      // Dynamic import for RDKit WASM
      const rdkitModule = await import('@rdkit/rdkit')
      const initRDKitModule = rdkitModule.default as (options?: object) => Promise<RDKitModule>
      rdkitInstance = await initRDKitModule()
      console.log('[RDKit] WASM loaded, version:', rdkitInstance.version())
      return rdkitInstance
    } catch (err) {
      console.error('[RDKit] Failed to load WASM:', err)
      throw err
    }
  })()

  return rdkitLoadPromise
}

interface MoleculeViewer2DProps {
  smiles: string
  title?: string
  highlights?: ScaffoldHighlight[]
  width?: number
  height?: number
  showControls?: boolean
  showSmiles?: boolean
  showProperties?: boolean
  onScaffoldClick?: (scaffold: ScaffoldHighlight) => void
  className?: string
}

export interface ScaffoldHighlight {
  id: string
  atomIndices: number[]
  color: string  // hex color
  label: string
  type: 'scaffold' | 'toxicity' | 'modification'
  smarts?: string  // SMARTS pattern for substructure matching
}

export interface MoleculeProperties {
  mw: number | null
  formula: string | null
  numAtoms: number | null
  numBonds: number | null
  numRings: number | null
  isValid: boolean
}

// SVG rendering options
const DEFAULT_SVG_OPTIONS = {
  width: 400,
  height: 300,
  bondLineWidth: 1.5,
  addAtomIndices: false,
  addStereoAnnotation: true,
  legendFontSize: 12,
  annotationFontScale: 0.7
}

export function MoleculeViewer2D({
  smiles,
  title = 'Molecular Structure',
  highlights = [],
  width = 400,
  height = 300,
  showControls = true,
  showSmiles = true,
  showProperties = false,
  onScaffoldClick,
  className
}: MoleculeViewer2DProps) {
  const containerRef = useRef<HTMLDivElement>(null)
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)
  const [copied, setCopied] = useState(false)
  const [scale, setScale] = useState(1)
  const [svgContent, setSvgContent] = useState<string>('')
  const [rdkit, setRdkit] = useState<RDKitModule | null>(null)
  const [properties, setProperties] = useState<MoleculeProperties | null>(null)
  const [activeHighlight, setActiveHighlight] = useState<string | null>(null)

  // Initialize RDKit
  useEffect(() => {
    initRDKit()
      .then(setRdkit)
      .catch(() => setError('Failed to load RDKit WASM module'))
  }, [])

  // Draw molecule with RDKit
  const drawMolecule = useCallback(async () => {
    if (!rdkit || !smiles) return

    setIsLoading(true)
    setError(null)

    let mol: JSMol | null = null

    try {
      mol = rdkit.get_mol(smiles)

      if (!mol || !mol.is_valid()) {
        setError('Invalid SMILES structure')
        setIsLoading(false)
        return
      }

      // Calculate scaled dimensions
      const scaledWidth = Math.round(width * scale)
      const scaledHeight = Math.round(height * scale)

      // Build highlight details for RDKit
      let svg: string

      if (highlights.length > 0 && activeHighlight) {
        const highlight = highlights.find(h => h.id === activeHighlight)
        if (highlight?.smarts) {
          // Use SMARTS-based substructure highlighting
          const qmol = rdkit.get_qmol(highlight.smarts)
          if (qmol && qmol.is_valid()) {
            const matchJson = mol.get_substruct_match(qmol)
            const match = JSON.parse(matchJson)

            if (match.atoms && match.atoms.length > 0) {
              const highlightDetails = JSON.stringify({
                width: scaledWidth,
                height: scaledHeight,
                atoms: match.atoms,
                bonds: match.bonds || [],
                highlightColour: hexToRgb(highlight.color),
                highlightRadius: 0.3
              })
              svg = mol.get_svg_with_highlights(highlightDetails)
            } else {
              svg = mol.get_svg(scaledWidth, scaledHeight)
            }
            qmol.delete()
          } else {
            svg = mol.get_svg(scaledWidth, scaledHeight)
          }
        } else if (highlight?.atomIndices.length) {
          // Use atom indices directly
          const highlightDetails = JSON.stringify({
            width: scaledWidth,
            height: scaledHeight,
            atoms: highlight.atomIndices,
            highlightColour: hexToRgb(highlight.color),
            highlightRadius: 0.3
          })
          svg = mol.get_svg_with_highlights(highlightDetails)
        } else {
          svg = mol.get_svg(scaledWidth, scaledHeight)
        }
      } else {
        svg = mol.get_svg(scaledWidth, scaledHeight)
      }

      // Clean up SVG for proper display
      svg = svg.replace(/width=['"]?\d+['"]?/, `width="${scaledWidth}"`)
      svg = svg.replace(/height=['"]?\d+['"]?/, `height="${scaledHeight}"`)

      setSvgContent(svg)

      // Extract properties if requested
      if (showProperties) {
        setProperties({
          mw: null, // RDKit.js basic doesn't expose MW directly
          formula: null,
          numAtoms: null,
          numBonds: null,
          numRings: null,
          isValid: true
        })
      }

      setIsLoading(false)
    } catch (err) {
      console.error('[RDKit] Draw error:', err)
      setError('Failed to render molecule')
      setIsLoading(false)
    } finally {
      if (mol) mol.delete()
    }
  }, [rdkit, smiles, width, height, scale, highlights, activeHighlight, showProperties])

  useEffect(() => {
    if (rdkit) {
      drawMolecule()
    }
  }, [drawMolecule, rdkit])

  const handleZoomIn = () => setScale(prev => Math.min(prev + 0.2, 2))
  const handleZoomOut = () => setScale(prev => Math.max(prev - 0.2, 0.5))
  const handleReset = () => {
    setScale(1)
    setActiveHighlight(null)
  }

  const handleDownload = () => {
    if (!svgContent) return

    // Create SVG blob and download
    const blob = new Blob([svgContent], { type: 'image/svg+xml' })
    const url = URL.createObjectURL(blob)
    const link = document.createElement('a')
    link.download = `molecule_${Date.now()}.svg`
    link.href = url
    link.click()
    URL.revokeObjectURL(url)
  }

  const handleDownloadPng = async () => {
    if (!svgContent || !containerRef.current) return

    // Convert SVG to PNG using canvas
    const canvas = document.createElement('canvas')
    const ctx = canvas.getContext('2d')
    if (!ctx) return

    const scaledWidth = Math.round(width * scale)
    const scaledHeight = Math.round(height * scale)
    canvas.width = scaledWidth
    canvas.height = scaledHeight

    const img = new Image()
    const svgBlob = new Blob([svgContent], { type: 'image/svg+xml;charset=utf-8' })
    const url = URL.createObjectURL(svgBlob)

    img.onload = () => {
      ctx.fillStyle = '#ffffff'
      ctx.fillRect(0, 0, scaledWidth, scaledHeight)
      ctx.drawImage(img, 0, 0)

      const pngUrl = canvas.toDataURL('image/png')
      const link = document.createElement('a')
      link.download = `molecule_${Date.now()}.png`
      link.href = pngUrl
      link.click()
      URL.revokeObjectURL(url)
    }
    img.src = url
  }

  const handleCopySmiles = async () => {
    await navigator.clipboard.writeText(smiles)
    setCopied(true)
    setTimeout(() => setCopied(false), 2000)
  }

  const handleHighlightClick = (highlight: ScaffoldHighlight) => {
    setActiveHighlight(prev => prev === highlight.id ? null : highlight.id)
    onScaffoldClick?.(highlight)
  }

  return (
    <Card className={className}>
      <CardHeader className="pb-2">
        <div className="flex items-center justify-between">
          <CardTitle className="text-sm font-medium flex items-center gap-2">
            <Atom className="w-4 h-4 text-purple-500" />
            {title}
            <Badge variant="secondary" className="text-[9px] px-1 py-0 h-4 bg-purple-100 text-purple-700">
              RDKit
            </Badge>
          </CardTitle>
          {showControls && (
            <div className="flex items-center gap-1">
              <Button
                variant="ghost"
                size="icon"
                className="h-7 w-7"
                onClick={handleZoomOut}
                title="Zoom Out"
              >
                <ZoomOut className="w-3 h-3" />
              </Button>
              <Badge variant="secondary" className="text-[10px] px-1">
                {(scale * 100).toFixed(0)}%
              </Badge>
              <Button
                variant="ghost"
                size="icon"
                className="h-7 w-7"
                onClick={handleZoomIn}
                title="Zoom In"
              >
                <ZoomIn className="w-3 h-3" />
              </Button>
              <Button
                variant="ghost"
                size="icon"
                className="h-7 w-7"
                onClick={handleReset}
                title="Reset"
              >
                <RotateCcw className="w-3 h-3" />
              </Button>
              <Button
                variant="ghost"
                size="icon"
                className="h-7 w-7"
                onClick={handleDownloadPng}
                title="Download PNG"
              >
                <Download className="w-3 h-3" />
              </Button>
            </div>
          )}
        </div>
      </CardHeader>
      <CardContent className="flex flex-col items-center">
        {/* SVG Container */}
        <div
          ref={containerRef}
          className="relative border rounded-lg bg-white overflow-hidden flex items-center justify-center"
          style={{
            width: width * scale,
            height: height * scale,
            minHeight: 100
          }}
        >
          {isLoading && (
            <div className="absolute inset-0 flex items-center justify-center bg-white/80 z-10">
              <div className="flex flex-col items-center gap-2">
                <Skeleton className="w-3/4 h-3/4" />
                <span className="text-xs text-gray-500">Loading RDKit...</span>
              </div>
            </div>
          )}
          {error && (
            <div className="absolute inset-0 flex flex-col items-center justify-center bg-red-50 z-10">
              <AlertCircle className="w-8 h-8 text-red-400 mb-2" />
              <p className="text-red-500 text-sm">{error}</p>
            </div>
          )}
          {svgContent && !error && (
            <div
              dangerouslySetInnerHTML={{ __html: svgContent }}
              className="molecule-svg"
              style={{ cursor: highlights.length > 0 ? 'pointer' : 'default' }}
            />
          )}
        </div>

        {/* Scaffold Legend / Highlight Controls */}
        {highlights.length > 0 && (
          <div className="flex flex-wrap gap-2 mt-2">
            {highlights.map((h) => (
              <Badge
                key={h.id}
                variant={activeHighlight === h.id ? 'default' : 'outline'}
                className="text-[10px] cursor-pointer hover:opacity-80 transition-all"
                style={{
                  borderColor: h.color,
                  color: activeHighlight === h.id ? '#fff' : h.color,
                  backgroundColor: activeHighlight === h.id ? h.color : 'transparent'
                }}
                onClick={() => handleHighlightClick(h)}
              >
                {h.label}
              </Badge>
            ))}
          </div>
        )}

        {/* Molecule Properties */}
        {showProperties && properties && (
          <div className="w-full mt-2 grid grid-cols-3 gap-2 text-[10px] text-gray-600">
            <div className="bg-slate-50 p-1 rounded text-center">
              <span className="font-medium">Valid:</span> {properties.isValid ? 'Yes' : 'No'}
            </div>
          </div>
        )}

        {/* SMILES Display */}
        {showSmiles && smiles && (
          <div className="w-full mt-3 pt-2 border-t">
            <div className="flex items-center justify-between gap-2">
              <code className="text-[10px] text-gray-600 font-mono break-all flex-1 bg-slate-50 p-2 rounded">
                {smiles.length > 80 ? `${smiles.slice(0, 80)}...` : smiles}
              </code>
              <Button
                variant="ghost"
                size="icon"
                className="h-6 w-6 shrink-0"
                onClick={handleCopySmiles}
                title="Copy SMILES"
              >
                {copied ? (
                  <Check className="w-3 h-3 text-green-500" />
                ) : (
                  <Copy className="w-3 h-3" />
                )}
              </Button>
            </div>
          </div>
        )}
      </CardContent>
    </Card>
  )
}

/**
 * Molecule Comparison View - 두 분자 나란히 비교
 */
interface MoleculeDiffViewerProps {
  originalSmiles: string
  optimizedSmiles: string
  originalLabel?: string
  optimizedLabel?: string
  toxicityHighlights?: ScaffoldHighlight[]
  improvementHighlights?: ScaffoldHighlight[]
}

export function MoleculeDiffViewer({
  originalSmiles,
  optimizedSmiles,
  originalLabel = 'Original',
  optimizedLabel = 'Optimized',
  toxicityHighlights = [],
  improvementHighlights = []
}: MoleculeDiffViewerProps) {
  return (
    <div className="grid grid-cols-2 gap-4">
      <div className="space-y-2">
        <Badge variant="outline" className="text-red-600 border-red-200">
          {originalLabel}
        </Badge>
        <MoleculeViewer2D
          smiles={originalSmiles}
          title=""
          highlights={toxicityHighlights}
          width={250}
          height={200}
          showControls={false}
          showSmiles={false}
        />
      </div>
      <div className="space-y-2">
        <Badge variant="outline" className="text-green-600 border-green-200">
          {optimizedLabel}
        </Badge>
        <MoleculeViewer2D
          smiles={optimizedSmiles}
          title=""
          highlights={improvementHighlights}
          width={250}
          height={200}
          showControls={false}
          showSmiles={false}
        />
      </div>
    </div>
  )
}

/**
 * Helper: Convert hex color to RGB array for RDKit
 */
function hexToRgb(hex: string): [number, number, number] {
  const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex)
  if (result) {
    return [
      parseInt(result[1], 16) / 255,
      parseInt(result[2], 16) / 255,
      parseInt(result[3], 16) / 255
    ]
  }
  return [1, 0.8, 0] // Default yellow
}

/**
 * Hook: Use RDKit for molecule operations
 */
export function useRDKit() {
  const [rdkit, setRdkit] = useState<RDKitModule | null>(null)
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    initRDKit()
      .then((module) => {
        setRdkit(module)
        setIsLoading(false)
      })
      .catch((err) => {
        setError(err.message)
        setIsLoading(false)
      })
  }, [])

  const validateSmiles = useCallback((smiles: string): boolean => {
    if (!rdkit) return false
    const mol = rdkit.get_mol(smiles)
    const isValid = mol ? mol.is_valid() : false
    if (mol) mol.delete()
    return isValid
  }, [rdkit])

  const canonicalizeSmiles = useCallback((smiles: string): string | null => {
    if (!rdkit) return null
    const mol = rdkit.get_mol(smiles)
    if (!mol || !mol.is_valid()) {
      if (mol) mol.delete()
      return null
    }
    const canonical = mol.get_smiles()
    mol.delete()
    return canonical
  }, [rdkit])

  const getMolBlock = useCallback((smiles: string): string | null => {
    if (!rdkit) return null
    const mol = rdkit.get_mol(smiles)
    if (!mol || !mol.is_valid()) {
      if (mol) mol.delete()
      return null
    }
    const molBlock = mol.get_molblock()
    mol.delete()
    return molBlock
  }, [rdkit])

  const hasSubstructMatch = useCallback((smiles: string, smarts: string): boolean => {
    if (!rdkit) return false
    const mol = rdkit.get_mol(smiles)
    const qmol = rdkit.get_qmol(smarts)

    if (!mol || !qmol || !mol.is_valid() || !qmol.is_valid()) {
      if (mol) mol.delete()
      if (qmol) qmol.delete()
      return false
    }

    const matchJson = mol.get_substruct_match(qmol)
    const match = JSON.parse(matchJson)

    mol.delete()
    qmol.delete()

    return match.atoms && match.atoms.length > 0
  }, [rdkit])

  return {
    rdkit,
    isLoading,
    error,
    validateSmiles,
    canonicalizeSmiles,
    getMolBlock,
    hasSubstructMatch
  }
}
