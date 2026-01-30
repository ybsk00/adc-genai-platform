/**
 * Design Runs Tab - Closed-loop R&D Pipeline UI
 * Displays design runs with frozen_params and assay results
 */
import { useState, useEffect } from 'react';
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table";
import { Dialog, DialogContent, DialogDescription, DialogHeader, DialogTitle, DialogTrigger } from "@/components/ui/dialog";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { RefreshCw, Eye, FlaskConical, CheckCircle2, XCircle, Clock, FileJson } from 'lucide-react';
import { toast } from "sonner";

const API_BASE = import.meta.env.VITE_API_BASE_URL || 'https://adc-backend-962229188169.asia-northeast3.run.app';

interface DesignRun {
    id: string;
    project_id: string;
    status: string;
    frozen_params: Record<string, unknown>;
    created_at: string;
    updated_at: string;
}

interface AssayResult {
    id: string;
    run_id: string;
    molecule_id: string;
    assay_type: string;
    raw_data: Record<string, unknown>;
    is_success: boolean;
    acceptance_criteria: Record<string, unknown>;
    confidence_score: number;
    created_at: string;
}

interface AssaySummary {
    total: number;
    success_count: number;
    failure_count: number;
    avg_confidence: number;
}

export function DesignRunsTab() {
    const [runs, setRuns] = useState<DesignRun[]>([]);
    const [selectedRun, setSelectedRun] = useState<DesignRun | null>(null);
    const [assayResults, setAssayResults] = useState<AssayResult[]>([]);
    const [assaySummary, setAssaySummary] = useState<AssaySummary | null>(null);
    const [loading, setLoading] = useState(false);

    const fetchRuns = async () => {
        setLoading(true);
        try {
            const res = await fetch(`${API_BASE}/api/design-runs`);
            const data = await res.json();
            setRuns(data);
        } catch (error) {
            console.error('Failed to fetch design runs:', error);
            toast.error("Failed to fetch design runs");
        } finally {
            setLoading(false);
        }
    };

    const fetchAssayResults = async (runId: string) => {
        try {
            const [resultsRes, summaryRes] = await Promise.all([
                fetch(`${API_BASE}/api/assay-results/by-run/${runId}`),
                fetch(`${API_BASE}/api/assay-results/summary/${runId}`)
            ]);
            const results = await resultsRes.json();
            const summary = await summaryRes.json();
            setAssayResults(results);
            setAssaySummary(summary);
        } catch (error) {
            console.error('Failed to fetch assay results:', error);
        }
    };

    useEffect(() => {
        fetchRuns();
    }, []);

    useEffect(() => {
        if (selectedRun) {
            fetchAssayResults(selectedRun.id);
        }
    }, [selectedRun]);

    const getStatusBadge = (status: string) => {
        const variants: Record<string, "default" | "secondary" | "destructive" | "outline"> = {
            draft: "secondary",
            computing: "outline",
            completed: "default",
            failed: "destructive"
        };
        return <Badge variant={variants[status] || "secondary"}>{status}</Badge>;
    };

    const formatDate = (dateStr: string) => {
        return new Date(dateStr).toLocaleString('ko-KR', {
            year: 'numeric', month: '2-digit', day: '2-digit',
            hour: '2-digit', minute: '2-digit'
        });
    };

    return (
        <div className="space-y-6">
            {/* Header */}
            <div className="flex items-center justify-between">
                <div>
                    <h2 className="text-2xl font-bold">Design Runs</h2>
                    <p className="text-muted-foreground">Closed-loop R&D Pipeline 실행 기록</p>
                </div>
                <Button onClick={fetchRuns} disabled={loading} variant="outline">
                    <RefreshCw className={`mr-2 h-4 w-4 ${loading ? 'animate-spin' : ''}`} />
                    Refresh
                </Button>
            </div>

            {/* Main Content */}
            <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
                {/* Design Runs List */}
                <Card className="lg:col-span-2">
                    <CardHeader>
                        <CardTitle className="flex items-center gap-2">
                            <FlaskConical className="h-5 w-5" />
                            Design Runs
                        </CardTitle>
                        <CardDescription>클릭하여 상세 정보 확인</CardDescription>
                    </CardHeader>
                    <CardContent>
                        <Table>
                            <TableHeader>
                                <TableRow>
                                    <TableHead>ID</TableHead>
                                    <TableHead>Status</TableHead>
                                    <TableHead>Version</TableHead>
                                    <TableHead>Created</TableHead>
                                    <TableHead>Actions</TableHead>
                                </TableRow>
                            </TableHeader>
                            <TableBody>
                                {runs.length === 0 ? (
                                    <TableRow>
                                        <TableCell colSpan={5} className="text-center text-muted-foreground py-8">
                                            No design runs found
                                        </TableCell>
                                    </TableRow>
                                ) : (
                                    runs.map((run) => (
                                        <TableRow
                                            key={run.id}
                                            className={`cursor-pointer ${selectedRun?.id === run.id ? 'bg-muted' : ''}`}
                                            onClick={() => setSelectedRun(run)}
                                        >
                                            <TableCell className="font-mono text-xs">{run.id.slice(0, 8)}...</TableCell>
                                            <TableCell>{getStatusBadge(run.status)}</TableCell>
                                            <TableCell>
                                                <Badge variant="outline">
                                                    {(run.frozen_params as { version?: string })?.version || 'v0.1'}
                                                </Badge>
                                            </TableCell>
                                            <TableCell className="text-sm">{formatDate(run.created_at)}</TableCell>
                                            <TableCell>
                                                <Dialog>
                                                    <DialogTrigger asChild>
                                                        <Button variant="ghost" size="sm">
                                                            <FileJson className="h-4 w-4" />
                                                        </Button>
                                                    </DialogTrigger>
                                                    <DialogContent className="max-w-2xl max-h-[80vh] overflow-auto">
                                                        <DialogHeader>
                                                            <DialogTitle>Frozen Parameters</DialogTitle>
                                                            <DialogDescription>
                                                                실행 당시의 스코링 파라미터 스냅샷 (Audit용)
                                                            </DialogDescription>
                                                        </DialogHeader>
                                                        <pre className="bg-muted p-4 rounded-lg text-xs overflow-auto">
                                                            {JSON.stringify(run.frozen_params, null, 2)}
                                                        </pre>
                                                    </DialogContent>
                                                </Dialog>
                                            </TableCell>
                                        </TableRow>
                                    ))
                                )}
                            </TableBody>
                        </Table>
                    </CardContent>
                </Card>

                {/* Selected Run Details */}
                <Card>
                    <CardHeader>
                        <CardTitle className="flex items-center gap-2">
                            <Eye className="h-5 w-5" />
                            Run Details
                        </CardTitle>
                    </CardHeader>
                    <CardContent>
                        {selectedRun ? (
                            <div className="space-y-4">
                                <div>
                                    <p className="text-sm text-muted-foreground">Run ID</p>
                                    <p className="font-mono text-xs">{selectedRun.id}</p>
                                </div>
                                <div>
                                    <p className="text-sm text-muted-foreground">Status</p>
                                    <p>{getStatusBadge(selectedRun.status)}</p>
                                </div>
                                <div>
                                    <p className="text-sm text-muted-foreground">Scoring Version</p>
                                    <p>{(selectedRun.frozen_params as { version?: string })?.version || 'N/A'}</p>
                                </div>
                                {assaySummary && (
                                    <div className="pt-4 border-t">
                                        <p className="text-sm font-medium mb-2">Assay Summary</p>
                                        <div className="grid grid-cols-2 gap-2 text-sm">
                                            <div className="flex items-center gap-1">
                                                <Clock className="h-3 w-3" />
                                                Total: {assaySummary.total}
                                            </div>
                                            <div className="flex items-center gap-1 text-green-600">
                                                <CheckCircle2 className="h-3 w-3" />
                                                Pass: {assaySummary.success_count}
                                            </div>
                                            <div className="flex items-center gap-1 text-red-600">
                                                <XCircle className="h-3 w-3" />
                                                Fail: {assaySummary.failure_count}
                                            </div>
                                            <div>
                                                Confidence: {(assaySummary.avg_confidence * 100).toFixed(1)}%
                                            </div>
                                        </div>
                                    </div>
                                )}
                            </div>
                        ) : (
                            <p className="text-muted-foreground text-center py-8">
                                Select a run to view details
                            </p>
                        )}
                    </CardContent>
                </Card>
            </div>

            {/* Assay Results Table */}
            {selectedRun && assayResults.length > 0 && (
                <Card>
                    <CardHeader>
                        <CardTitle>Assay Results</CardTitle>
                        <CardDescription>Run ID: {selectedRun.id.slice(0, 8)}...</CardDescription>
                    </CardHeader>
                    <CardContent>
                        <Table>
                            <TableHeader>
                                <TableRow>
                                    <TableHead>Assay Type</TableHead>
                                    <TableHead>Molecule</TableHead>
                                    <TableHead>Result</TableHead>
                                    <TableHead>Criteria</TableHead>
                                    <TableHead>Confidence</TableHead>
                                    <TableHead>Date</TableHead>
                                </TableRow>
                            </TableHeader>
                            <TableBody>
                                {assayResults.map((result) => (
                                    <TableRow key={result.id}>
                                        <TableCell>
                                            <Badge variant="outline">{result.assay_type}</Badge>
                                        </TableCell>
                                        <TableCell className="font-mono text-xs">
                                            {result.molecule_id.slice(0, 8)}...
                                        </TableCell>
                                        <TableCell>
                                            {result.is_success ? (
                                                <Badge className="bg-green-500">Pass</Badge>
                                            ) : (
                                                <Badge variant="destructive">Fail</Badge>
                                            )}
                                        </TableCell>
                                        <TableCell className="text-xs">
                                            {JSON.stringify(result.acceptance_criteria).slice(0, 30)}...
                                        </TableCell>
                                        <TableCell>
                                            <span className={result.confidence_score > 0.7 ? 'text-green-600' : 'text-yellow-600'}>
                                                {(result.confidence_score * 100).toFixed(1)}%
                                            </span>
                                        </TableCell>
                                        <TableCell className="text-sm">{formatDate(result.created_at)}</TableCell>
                                    </TableRow>
                                ))}
                            </TableBody>
                        </Table>
                    </CardContent>
                </Card>
            )}
        </div>
    );
}
