/**
 * MolStar 3D Viewer Component
 * PDB/분자 구조를 3D로 시각화
 * 
 * [Dev Note]
 * Mol* (MolStar)는 WebGL 기반 분자 시각화 라이브러리
 * PDB ID 또는 URL을 입력받아 3D 구조를 렌더링
 */
import { useEffect, useRef, useState } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import {
    RotateCcw,
    ZoomIn,
    ZoomOut,
    Maximize2,
    Download,
    Loader2
} from 'lucide-react'

interface MolStarViewerProps {
    pdbId?: string        // PDB ID (e.g., "1HZH")
    pdbUrl?: string       // Direct PDB file URL
    height?: number       // Height in pixels
    showControls?: boolean
    onLoad?: () => void
    onError?: (error: string) => void
}

export function MolStarViewer({
    pdbId,
    pdbUrl,
    height = 400,
    showControls = true,
    onLoad,
    onError
}: MolStarViewerProps) {
    const containerRef = useRef<HTMLDivElement>(null)
    const [isLoading, setIsLoading] = useState(true)
    const [error, setError] = useState<string | null>(null)

    // Mol* plugin instance (stored in ref to persist across renders)
    const pluginRef = useRef<any>(null)

    useEffect(() => {
        if (!containerRef.current) return

        // Load Mol* dynamically
        const loadMolStar = async () => {
            setIsLoading(true)
            setError(null)

            try {
                // Dynamic import of Mol* library
                // @ts-ignore
                const { createPluginUI } = await import('molstar/lib/mol-plugin-ui')
                // @ts-ignore
                const { DefaultPluginUISpec } = await import('molstar/lib/mol-plugin-ui/spec')

                // Create plugin instance
                const plugin = await createPluginUI(containerRef.current!, {
                    ...DefaultPluginUISpec(),
                    layout: {
                        initial: {
                            isExpanded: false,
                            showControls: false,
                            controlsDisplay: 'reactive'
                        }
                    }
                })

                pluginRef.current = plugin

                // Load structure
                const url = pdbUrl || `https://files.rcsb.org/download/${pdbId}.pdb`

                await plugin.dataTransaction(async () => {
                    const data = await plugin.builders.data.download({
                        url,
                        isBinary: false
                    })

                    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb')
                    await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default')
                })

                setIsLoading(false)
                onLoad?.()

            } catch (err) {
                const errorMessage = err instanceof Error ? err.message : 'Failed to load structure'
                setError(errorMessage)
                setIsLoading(false)
                onError?.(errorMessage)
            }
        }

        loadMolStar()

        return () => {
            // Cleanup
            if (pluginRef.current) {
                pluginRef.current.dispose()
            }
        }
    }, [pdbId, pdbUrl])

    const handleReset = () => {
        pluginRef.current?.canvas3d?.requestCameraReset()
    }

    const handleZoomIn = () => {
        // Zoom in by adjusting camera distance
        const camera = pluginRef.current?.canvas3d?.camera
        if (camera) {
            // Implementation depends on Mol* version
        }
    }

    const handleZoomOut = () => {
        // Zoom out
    }

    const handleDownloadImage = () => {
        pluginRef.current?.helpers.viewportScreenshot?.download()
    }

    // Fallback: Show placeholder with mock 3D preview
    const renderPlaceholder = () => (
        <div
            className="relative bg-gradient-to-br from-slate-900 to-slate-800 rounded-lg overflow-hidden"
            style={{ height }}
        >
            {/* Mock 3D structure preview */}
            <div className="absolute inset-0 flex items-center justify-center">
                <svg viewBox="0 0 200 200" className="w-48 h-48 opacity-50">
                    <defs>
                        <linearGradient id="molGrad" x1="0%" y1="0%" x2="100%" y2="100%">
                            <stop offset="0%" stopColor="#6366f1" />
                            <stop offset="100%" stopColor="#ec4899" />
                        </linearGradient>
                    </defs>
                    {/* Simplified ADC structure */}
                    <g fill="none" stroke="url(#molGrad)" strokeWidth="2">
                        {/* Y-shaped antibody */}
                        <path d="M100 180 L100 120 M100 120 L60 80 M100 120 L140 80" strokeWidth="4" />
                        <circle cx="60" cy="80" r="15" fill="url(#molGrad)" fillOpacity="0.3" />
                        <circle cx="140" cy="80" r="15" fill="url(#molGrad)" fillOpacity="0.3" />
                        <circle cx="100" cy="120" r="10" fill="url(#molGrad)" fillOpacity="0.5" />

                        {/* Payload dots */}
                        <circle cx="45" cy="65" r="8" fill="#ec4899" />
                        <circle cx="155" cy="65" r="8" fill="#ec4899" />
                        <circle cx="60" cy="95" r="6" fill="#ec4899" />
                        <circle cx="140" cy="95" r="6" fill="#ec4899" />

                        {/* Linker lines */}
                        <line x1="52" y1="72" x2="60" y2="80" stroke="#ec4899" strokeWidth="1.5" />
                        <line x1="148" y1="72" x2="140" y2="80" stroke="#ec4899" strokeWidth="1.5" />
                    </g>
                </svg>
            </div>

            {/* Loading overlay */}
            {isLoading && (
                <div className="absolute inset-0 bg-slate-900/80 flex items-center justify-center">
                    <div className="text-center">
                        <Loader2 className="w-8 h-8 text-purple-500 animate-spin mx-auto" />
                        <p className="text-sm text-slate-400 mt-2">구조 로딩 중...</p>
                    </div>
                </div>
            )}

            {/* Info badge */}
            <div className="absolute bottom-3 left-3 bg-slate-800/80 backdrop-blur-sm rounded-lg px-3 py-1.5">
                <p className="text-xs text-slate-300">
                    {pdbId ? `PDB: ${pdbId}` : 'ADC Structure Preview'}
                </p>
            </div>
        </div>
    )

    return (
        <div className="space-y-2">
            {/* Viewer */}
            <div
                ref={containerRef}
                className="relative rounded-lg overflow-hidden"
                style={{ height }}
            >
                {/* If Mol* fails to load, show placeholder */}
                {(error || !pdbId && !pdbUrl) && renderPlaceholder()}
            </div>

            {/* Controls */}
            {showControls && (
                <div className="flex items-center justify-between">
                    <div className="flex gap-1">
                        <Button
                            variant="outline"
                            size="sm"
                            onClick={handleReset}
                            className="h-8"
                        >
                            <RotateCcw className="w-4 h-4" />
                        </Button>
                        <Button
                            variant="outline"
                            size="sm"
                            onClick={handleZoomIn}
                            className="h-8"
                        >
                            <ZoomIn className="w-4 h-4" />
                        </Button>
                        <Button
                            variant="outline"
                            size="sm"
                            onClick={handleZoomOut}
                            className="h-8"
                        >
                            <ZoomOut className="w-4 h-4" />
                        </Button>
                    </div>
                    <div className="flex gap-1">
                        <Button
                            variant="outline"
                            size="sm"
                            onClick={handleDownloadImage}
                            className="h-8"
                        >
                            <Download className="w-4 h-4 mr-1" />
                            이미지
                        </Button>
                        <Button
                            variant="outline"
                            size="sm"
                            className="h-8"
                        >
                            <Maximize2 className="w-4 h-4" />
                        </Button>
                    </div>
                </div>
            )}
        </div>
    )
}

/**
 * Standalone 3D Viewer Card
 * Dashboard에서 독립적으로 사용 가능
 */
export function MolStarViewerCard({
    pdbId = "1HZH",
    title = "3D Structure"
}: {
    pdbId?: string
    title?: string
}) {
    return (
        <Card>
            <CardHeader className="pb-2">
                <CardTitle className="text-lg">{title}</CardTitle>
            </CardHeader>
            <CardContent>
                <MolStarViewer pdbId={pdbId} height={350} />
            </CardContent>
        </Card>
    )
}
