import { Sheet, SheetContent, SheetHeader, SheetTitle, SheetDescription } from '@/components/ui/sheet'
import { Badge } from '@/components/ui/badge'
import { Button } from '@/components/ui/button'
import { FileDown, ExternalLink } from 'lucide-react'
import { format } from 'date-fns'

interface SimulationLog {
    id: string
    title: string
    status: 'completed' | 'failed' | 'processing' | 'pending'
    created_at: string
    report_url: string | null
    input_params: {
        target?: string
        antibody?: string
        payload?: string
        linker?: string
        dar?: number
    } | null
    result_summary: {
        success_score?: number
        toxicity_warnings?: string[]
        binding_affinity?: number
        dar_prediction?: number
        grade?: string
        recommendation?: string
    } | null
    profiles: {
        email: string
    } | null
}

interface SimulationDetailDrawerProps {
    simulation: SimulationLog | null
    open: boolean
    onClose: () => void
}

export function SimulationDetailDrawer({ simulation, open, onClose }: SimulationDetailDrawerProps) {
    if (!simulation) return null

    const handleDownloadPDF = () => {
        if (simulation.report_url) {
            window.open(simulation.report_url, '_blank')
        }
    }

    return (
        <Sheet open={open} onOpenChange={onClose}>
            <SheetContent className="bg-slate-900 border-slate-800 text-white w-full sm:max-w-2xl overflow-y-auto">
                <SheetHeader>
                    <SheetTitle className="text-white">{simulation.title}</SheetTitle>
                    <SheetDescription className="text-slate-400">
                        Simulation Details • {format(new Date(simulation.created_at), 'PPpp')}
                    </SheetDescription>
                </SheetHeader>

                <div className="mt-6 space-y-6">
                    {/* Section A: Input Parameters */}
                    <div className="space-y-3">
                        <h3 className="text-lg font-semibold text-white border-b border-slate-700 pb-2">
                            Input Parameters
                        </h3>
                        <div className="grid grid-cols-2 gap-4">
                            <div>
                                <p className="text-sm text-slate-400">Target Protein</p>
                                <p className="text-white font-medium">
                                    {simulation.input_params?.target || 'N/A'}
                                </p>
                            </div>
                            <div>
                                <p className="text-sm text-slate-400">Antibody</p>
                                <p className="text-white font-medium">
                                    {simulation.input_params?.antibody || 'N/A'}
                                </p>
                            </div>
                            <div>
                                <p className="text-sm text-slate-400">Payload</p>
                                <p className="text-white font-medium">
                                    {simulation.input_params?.payload || 'N/A'}
                                </p>
                            </div>
                            <div>
                                <p className="text-sm text-slate-400">Linker</p>
                                <p className="text-white font-medium">
                                    {simulation.input_params?.linker || 'N/A'}
                                </p>
                            </div>
                            {simulation.input_params?.dar && (
                                <div>
                                    <p className="text-sm text-slate-400">DAR (Drug-Antibody Ratio)</p>
                                    <p className="text-white font-medium">
                                        {simulation.input_params.dar}
                                    </p>
                                </div>
                            )}
                        </div>
                    </div>

                    {/* Section B: Simulation Results */}
                    {simulation.result_summary && (
                        <div className="space-y-3">
                            <h3 className="text-lg font-semibold text-white border-b border-slate-700 pb-2">
                                Simulation Results
                            </h3>
                            <div className="space-y-3">
                                {simulation.result_summary.success_score !== undefined && (
                                    <div>
                                        <p className="text-sm text-slate-400">Success Score</p>
                                        <div className="flex items-center gap-2">
                                            <div className="flex-1 bg-slate-800 rounded-full h-2">
                                                <div
                                                    className="bg-green-500 h-2 rounded-full"
                                                    style={{ width: `${simulation.result_summary.success_score}%` }}
                                                />
                                            </div>
                                            <span className="text-white font-medium">
                                                {simulation.result_summary.success_score}%
                                            </span>
                                        </div>
                                    </div>
                                )}

                                {simulation.result_summary.grade && (
                                    <div>
                                        <p className="text-sm text-slate-400">Grade</p>
                                        <Badge
                                            variant="outline"
                                            className={
                                                simulation.result_summary.grade === 'A'
                                                    ? 'border-green-500/30 text-green-400'
                                                    : simulation.result_summary.grade === 'B'
                                                        ? 'border-blue-500/30 text-blue-400'
                                                        : 'border-yellow-500/30 text-yellow-400'
                                            }
                                        >
                                            {simulation.result_summary.grade}
                                        </Badge>
                                    </div>
                                )}

                                {simulation.result_summary.recommendation && (
                                    <div>
                                        <p className="text-sm text-slate-400">Recommendation</p>
                                        <p className="text-white font-medium">
                                            {simulation.result_summary.recommendation}
                                        </p>
                                    </div>
                                )}

                                {simulation.result_summary.toxicity_warnings && simulation.result_summary.toxicity_warnings.length > 0 && (
                                    <div>
                                        <p className="text-sm text-slate-400 mb-2">Toxicity Warnings</p>
                                        <div className="space-y-1">
                                            {simulation.result_summary.toxicity_warnings.map((warning, idx) => (
                                                <div key={idx} className="flex items-start gap-2 text-sm">
                                                    <span className="text-red-400">⚠️</span>
                                                    <span className="text-slate-300">{warning}</span>
                                                </div>
                                            ))}
                                        </div>
                                    </div>
                                )}

                                {simulation.result_summary.binding_affinity !== undefined && (
                                    <div>
                                        <p className="text-sm text-slate-400">Binding Affinity (Kd)</p>
                                        <p className="text-white font-medium">
                                            {simulation.result_summary.binding_affinity} nM
                                        </p>
                                    </div>
                                )}

                                {simulation.result_summary.dar_prediction !== undefined && (
                                    <div>
                                        <p className="text-sm text-slate-400">DAR Prediction</p>
                                        <p className="text-white font-medium">
                                            {simulation.result_summary.dar_prediction}
                                        </p>
                                    </div>
                                )}
                            </div>
                        </div>
                    )}

                    {/* Section C: Artifacts */}
                    <div className="space-y-3">
                        <h3 className="text-lg font-semibold text-white border-b border-slate-700 pb-2">
                            Generated Report
                        </h3>
                        {simulation.report_url ? (
                            <div className="flex gap-2">
                                <Button
                                    onClick={handleDownloadPDF}
                                    className="bg-purple-500 hover:bg-purple-600 flex-1"
                                >
                                    <FileDown className="w-4 h-4 mr-2" />
                                    Download PDF Report
                                </Button>
                                <Button
                                    onClick={handleDownloadPDF}
                                    variant="outline"
                                    className="border-slate-700"
                                >
                                    <ExternalLink className="w-4 h-4" />
                                </Button>
                            </div>
                        ) : (
                            <p className="text-slate-500 text-sm">No report available</p>
                        )}
                    </div>

                    {/* User Info */}
                    <div className="pt-4 border-t border-slate-800">
                        <p className="text-sm text-slate-400">
                            User: <span className="text-white">{simulation.profiles?.email || 'Unknown'}</span>
                        </p>
                        <p className="text-sm text-slate-400">
                            Status: <Badge
                                variant="outline"
                                className={
                                    simulation.status === 'completed'
                                        ? 'border-green-500/30 text-green-400'
                                        : simulation.status === 'failed'
                                            ? 'border-red-500/30 text-red-400'
                                            : simulation.status === 'processing'
                                                ? 'border-blue-500/30 text-blue-400'
                                                : 'border-slate-600 text-slate-400'
                                }
                            >
                                {simulation.status}
                            </Badge>
                        </p>
                    </div>
                </div>
            </SheetContent>
        </Sheet>
    )
}
