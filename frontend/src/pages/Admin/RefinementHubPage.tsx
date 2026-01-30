/**
 * Data Refinement Hub Page
 * v2.2 Management - Batch AI Fixer & Data Quality Control
 */

import { useState, useEffect } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import {
  Wrench,
  RefreshCw,
  Sparkles,
  CheckCircle,
  XCircle,
  AlertTriangle,
  Database,
  FileSearch,
  Beaker,
  FlaskConical,
  Play,
  Pause,
  RotateCcw,
  Eye,
  Check,
  X,
  ChevronDown,
  ChevronRight,
  Filter
} from 'lucide-react';
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Progress } from '@/components/ui/progress';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { ScrollArea } from '@/components/ui/scroll-area';
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select';
import { toast } from 'sonner';
import { cn } from '@/lib/utils';
import { API_BASE_URL } from '@/lib/api';

// ============================================================================
// Types
// ============================================================================

interface QuarantinedItem {
  id: string;
  sourceTable: string;
  sourceId: string;
  errorType: string;
  errorDetails: string;
  originalData: Record<string, unknown>;
  suggestedFix: Record<string, unknown> | null;
  status: 'pending' | 'reviewing' | 'approved' | 'rejected';
  priority: 'low' | 'normal' | 'high' | 'critical';
  createdAt: string;
}

interface AIFixSuggestion {
  id: string;
  quarantinedId: string;
  suggestionType: string;
  originalValue: string;
  suggestedValue: string;
  aiConfidence: number;
  aiReasoning: string;
  status: 'pending' | 'applied' | 'rejected';
}

interface DataHealthStats {
  sourceTable: string;
  totalCount: number;
  normalizedCount: number;
  embeddedCount: number;
  normalizationRate: number;
  embeddingRate: number;
}

interface BatchJob {
  id: string;
  type: 'smiles_validation' | 'target_normalization' | 'embedding_generation';
  status: 'pending' | 'running' | 'completed' | 'failed';
  totalItems: number;
  processedItems: number;
  successCount: number;
  errorCount: number;
  startedAt: string | null;
  completedAt: string | null;
}

// ============================================================================
// Component
// ============================================================================

export function RefinementHubPage() {
  // State
  const [quarantinedItems, setQuarantinedItems] = useState<QuarantinedItem[]>([]);
  const [aiSuggestions, setAiSuggestions] = useState<AIFixSuggestion[]>([]);
  const [dataHealth, setDataHealth] = useState<DataHealthStats[]>([]);
  const [batchJobs, setBatchJobs] = useState<BatchJob[]>([]);
  const [selectedItem, setSelectedItem] = useState<QuarantinedItem | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [filterStatus, setFilterStatus] = useState<string>('all');
  const [filterTable, setFilterTable] = useState<string>('all');

  // Load data on mount
  useEffect(() => {
    loadQuarantinedData();
    loadDataHealth();
    loadBatchJobs();
  }, []);

  const loadQuarantinedData = async () => {
    try {
      const response = await fetch(`${API_BASE_URL}/api/admin/refinement/quarantined`);
      if (response.ok) {
        const data = await response.json();
        setQuarantinedItems(data.items || []);
      }
    } catch (error) {
      console.error('Failed to load quarantined data:', error);
      // Mock data for demo
      setQuarantinedItems([
        {
          id: '1',
          sourceTable: 'commercial_reagents',
          sourceId: 'uuid-1',
          errorType: 'invalid_smiles',
          errorDetails: 'SMILES string contains invalid characters',
          originalData: { smiles: 'CC(=O)Oc1ccccc1C(=O)O[X]', name: 'Aspirin derivative' },
          suggestedFix: { smiles: 'CC(=O)Oc1ccccc1C(=O)O' },
          status: 'pending',
          priority: 'high',
          createdAt: new Date().toISOString(),
        },
        {
          id: '2',
          sourceTable: 'antibody_library',
          sourceId: 'uuid-2',
          errorType: 'unknown_target',
          errorDetails: 'Target "HER-2/NEU" not in canonical list',
          originalData: { target: 'HER-2/NEU', name: 'Trastuzumab' },
          suggestedFix: { target: 'HER2', target_normalized: 'HER2' },
          status: 'pending',
          priority: 'normal',
          createdAt: new Date().toISOString(),
        },
        {
          id: '3',
          sourceTable: 'commercial_reagents',
          sourceId: 'uuid-3',
          errorType: 'missing_smiles',
          errorDetails: 'No SMILES data available for linker',
          originalData: { name: 'MC-VC-PAB', category: 'linker' },
          suggestedFix: null,
          status: 'pending',
          priority: 'critical',
          createdAt: new Date().toISOString(),
        },
      ]);
    }
  };

  const loadDataHealth = async () => {
    try {
      const response = await fetch(`${API_BASE_URL}/api/admin/refinement/health`);
      if (response.ok) {
        const data = await response.json();
        setDataHealth(data.metrics || []);
      }
    } catch (error) {
      console.error('Failed to load data health:', error);
      // Mock data
      setDataHealth([
        {
          sourceTable: 'antibody_library',
          totalCount: 4520,
          normalizedCount: 4350,
          embeddedCount: 3800,
          normalizationRate: 96.2,
          embeddingRate: 84.1,
        },
        {
          sourceTable: 'commercial_reagents',
          totalCount: 12500,
          normalizedCount: 11200,
          embeddedCount: 8500,
          normalizationRate: 89.6,
          embeddingRate: 68.0,
        },
        {
          sourceTable: 'golden_set',
          totalCount: 350,
          normalizedCount: 350,
          embeddedCount: 350,
          normalizationRate: 100,
          embeddingRate: 100,
        },
      ]);
    }
  };

  const loadBatchJobs = async () => {
    try {
      const response = await fetch(`${API_BASE_URL}/api/admin/refinement/jobs`);
      if (response.ok) {
        const data = await response.json();
        setBatchJobs(data.jobs || []);
      }
    } catch (error) {
      console.error('Failed to load batch jobs:', error);
      // Mock data
      setBatchJobs([
        {
          id: 'job-1',
          type: 'target_normalization',
          status: 'running',
          totalItems: 500,
          processedItems: 234,
          successCount: 230,
          errorCount: 4,
          startedAt: new Date(Date.now() - 300000).toISOString(),
          completedAt: null,
        },
        {
          id: 'job-2',
          type: 'smiles_validation',
          status: 'completed',
          totalItems: 1200,
          processedItems: 1200,
          successCount: 1150,
          errorCount: 50,
          startedAt: new Date(Date.now() - 3600000).toISOString(),
          completedAt: new Date(Date.now() - 1800000).toISOString(),
        },
      ]);
    }
  };

  // Handlers
  const handleRunBatchAIFix = async (type: string) => {
    setIsLoading(true);
    try {
      const response = await fetch(`${API_BASE_URL}/api/admin/refinement/batch-fix`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ type }),
      });

      if (response.ok) {
        toast.success('Batch AI fix job started');
        loadBatchJobs();
      } else {
        toast.error('Failed to start batch fix');
      }
    } catch (error) {
      toast.error('Failed to start batch fix');
    } finally {
      setIsLoading(false);
    }
  };

  const handleApproveItem = async (itemId: string) => {
    try {
      const response = await fetch(`${API_BASE_URL}/api/admin/refinement/quarantined/${itemId}/approve`, {
        method: 'POST',
      });

      if (response.ok) {
        toast.success('Fix approved and applied');
        setQuarantinedItems(prev => prev.filter(item => item.id !== itemId));
        setSelectedItem(null);
      }
    } catch (error) {
      toast.error('Failed to approve fix');
    }
  };

  const handleRejectItem = async (itemId: string) => {
    try {
      const response = await fetch(`${API_BASE_URL}/api/admin/refinement/quarantined/${itemId}/reject`, {
        method: 'POST',
      });

      if (response.ok) {
        toast.success('Item rejected');
        setQuarantinedItems(prev =>
          prev.map(item => item.id === itemId ? { ...item, status: 'rejected' as const } : item)
        );
      }
    } catch (error) {
      toast.error('Failed to reject item');
    }
  };

  const handleGenerateAISuggestion = async (itemId: string) => {
    setIsLoading(true);
    try {
      const response = await fetch(`${API_BASE_URL}/api/admin/refinement/quarantined/${itemId}/ai-suggest`, {
        method: 'POST',
      });

      if (response.ok) {
        const data = await response.json();
        toast.success('AI suggestion generated');
        loadQuarantinedData();
      }
    } catch (error) {
      toast.error('Failed to generate AI suggestion');
    } finally {
      setIsLoading(false);
    }
  };

  // Filtered items
  const filteredItems = quarantinedItems.filter(item => {
    if (filterStatus !== 'all' && item.status !== filterStatus) return false;
    if (filterTable !== 'all' && item.sourceTable !== filterTable) return false;
    return true;
  });

  // Stats
  const pendingCount = quarantinedItems.filter(i => i.status === 'pending').length;
  const criticalCount = quarantinedItems.filter(i => i.priority === 'critical').length;

  const getPriorityColor = (priority: string) => {
    switch (priority) {
      case 'critical': return 'bg-red-500/20 text-red-400 border-red-500/30';
      case 'high': return 'bg-amber-500/20 text-amber-400 border-amber-500/30';
      case 'normal': return 'bg-blue-500/20 text-blue-400 border-blue-500/30';
      default: return 'bg-slate-500/20 text-slate-400 border-slate-500/30';
    }
  };

  const getStatusColor = (status: string) => {
    switch (status) {
      case 'pending': return 'bg-amber-500/20 text-amber-400';
      case 'reviewing': return 'bg-blue-500/20 text-blue-400';
      case 'approved': return 'bg-emerald-500/20 text-emerald-400';
      case 'rejected': return 'bg-red-500/20 text-red-400';
      default: return 'bg-slate-500/20 text-slate-400';
    }
  };

  const getHealthColor = (rate: number) => {
    if (rate >= 95) return 'text-emerald-400';
    if (rate >= 80) return 'text-amber-400';
    return 'text-red-400';
  };

  return (
    <div className="space-y-6">
      {/* Page Header */}
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-2xl font-bold text-white flex items-center gap-3">
            <Wrench className="h-7 w-7 text-purple-400" />
            Data Refinement Hub
          </h1>
          <p className="text-slate-400 mt-1">
            AI-powered data quality control and batch corrections
          </p>
        </div>
        <div className="flex gap-2">
          <Button
            variant="outline"
            onClick={() => { loadQuarantinedData(); loadDataHealth(); loadBatchJobs(); }}
            className="border-slate-700"
          >
            <RefreshCw className="h-4 w-4 mr-2" />
            Refresh
          </Button>
        </div>
      </div>

      {/* Quick Stats */}
      <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
        <Card className="bg-slate-900/50 border-slate-700">
          <CardContent className="p-4">
            <div className="flex items-center gap-3">
              <div className="p-2 bg-amber-500/20 rounded-lg">
                <AlertTriangle className="h-5 w-5 text-amber-400" />
              </div>
              <div>
                <p className="text-2xl font-bold text-white">{pendingCount}</p>
                <p className="text-sm text-slate-400">Pending Review</p>
              </div>
            </div>
          </CardContent>
        </Card>

        <Card className="bg-slate-900/50 border-slate-700">
          <CardContent className="p-4">
            <div className="flex items-center gap-3">
              <div className="p-2 bg-red-500/20 rounded-lg">
                <XCircle className="h-5 w-5 text-red-400" />
              </div>
              <div>
                <p className="text-2xl font-bold text-white">{criticalCount}</p>
                <p className="text-sm text-slate-400">Critical Issues</p>
              </div>
            </div>
          </CardContent>
        </Card>

        <Card className="bg-slate-900/50 border-slate-700">
          <CardContent className="p-4">
            <div className="flex items-center gap-3">
              <div className="p-2 bg-emerald-500/20 rounded-lg">
                <Sparkles className="h-5 w-5 text-emerald-400" />
              </div>
              <div>
                <p className="text-2xl font-bold text-white">
                  {batchJobs.filter(j => j.status === 'running').length}
                </p>
                <p className="text-sm text-slate-400">Active Jobs</p>
              </div>
            </div>
          </CardContent>
        </Card>

        <Card className="bg-slate-900/50 border-slate-700">
          <CardContent className="p-4">
            <div className="flex items-center gap-3">
              <div className="p-2 bg-purple-500/20 rounded-lg">
                <Database className="h-5 w-5 text-purple-400" />
              </div>
              <div>
                <p className="text-2xl font-bold text-white">
                  {dataHealth.reduce((acc, d) => acc + d.totalCount, 0).toLocaleString()}
                </p>
                <p className="text-sm text-slate-400">Total Records</p>
              </div>
            </div>
          </CardContent>
        </Card>
      </div>

      <Tabs defaultValue="quarantine" className="space-y-6">
        <TabsList className="bg-slate-800">
          <TabsTrigger value="quarantine">Quarantine Queue</TabsTrigger>
          <TabsTrigger value="batch">Batch AI Fixer</TabsTrigger>
          <TabsTrigger value="health">Data Health</TabsTrigger>
        </TabsList>

        {/* Quarantine Queue Tab */}
        <TabsContent value="quarantine" className="space-y-6">
          <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
            {/* List */}
            <div className="lg:col-span-2">
              <Card className="bg-slate-900/50 border-slate-700">
                <CardHeader className="pb-3">
                  <div className="flex items-center justify-between">
                    <CardTitle className="text-white">Quarantined Items</CardTitle>
                    <div className="flex gap-2">
                      <Select value={filterStatus} onValueChange={setFilterStatus}>
                        <SelectTrigger className="w-[140px] bg-slate-800 border-slate-600">
                          <SelectValue placeholder="Status" />
                        </SelectTrigger>
                        <SelectContent>
                          <SelectItem value="all">All Status</SelectItem>
                          <SelectItem value="pending">Pending</SelectItem>
                          <SelectItem value="reviewing">Reviewing</SelectItem>
                          <SelectItem value="approved">Approved</SelectItem>
                          <SelectItem value="rejected">Rejected</SelectItem>
                        </SelectContent>
                      </Select>
                      <Select value={filterTable} onValueChange={setFilterTable}>
                        <SelectTrigger className="w-[180px] bg-slate-800 border-slate-600">
                          <SelectValue placeholder="Source" />
                        </SelectTrigger>
                        <SelectContent>
                          <SelectItem value="all">All Sources</SelectItem>
                          <SelectItem value="antibody_library">Antibody Library</SelectItem>
                          <SelectItem value="commercial_reagents">Commercial Reagents</SelectItem>
                          <SelectItem value="golden_set">Golden Set</SelectItem>
                        </SelectContent>
                      </Select>
                    </div>
                  </div>
                </CardHeader>
                <CardContent>
                  <ScrollArea className="h-[500px]">
                    <div className="space-y-2">
                      {filteredItems.length === 0 ? (
                        <div className="text-center py-12 text-slate-500">
                          <CheckCircle className="h-12 w-12 mx-auto mb-4 opacity-50" />
                          <p>No quarantined items</p>
                        </div>
                      ) : (
                        filteredItems.map((item) => (
                          <div
                            key={item.id}
                            onClick={() => setSelectedItem(item)}
                            className={cn(
                              "p-4 rounded-lg border cursor-pointer transition-all",
                              selectedItem?.id === item.id
                                ? "bg-purple-900/30 border-purple-500"
                                : "bg-slate-800/50 border-slate-700 hover:border-slate-600"
                            )}
                          >
                            <div className="flex items-start justify-between">
                              <div className="flex-1">
                                <div className="flex items-center gap-2 mb-1">
                                  <Badge variant="outline" className={getPriorityColor(item.priority)}>
                                    {item.priority}
                                  </Badge>
                                  <Badge className={getStatusColor(item.status)}>
                                    {item.status}
                                  </Badge>
                                </div>
                                <p className="text-white font-medium">{item.errorType.replace(/_/g, ' ')}</p>
                                <p className="text-sm text-slate-400 mt-1">{item.errorDetails}</p>
                                <p className="text-xs text-slate-500 mt-2">
                                  Source: {item.sourceTable}
                                </p>
                              </div>
                              <ChevronRight className={cn(
                                "h-5 w-5 text-slate-500 transition-transform",
                                selectedItem?.id === item.id && "rotate-90"
                              )} />
                            </div>
                          </div>
                        ))
                      )}
                    </div>
                  </ScrollArea>
                </CardContent>
              </Card>
            </div>

            {/* Detail Panel */}
            <div>
              <Card className="bg-slate-900/50 border-slate-700 sticky top-4">
                <CardHeader>
                  <CardTitle className="text-white">Item Details</CardTitle>
                </CardHeader>
                <CardContent>
                  {selectedItem ? (
                    <div className="space-y-4">
                      <div>
                        <Label className="text-slate-400 text-xs">Error Type</Label>
                        <p className="text-white font-medium">
                          {selectedItem.errorType.replace(/_/g, ' ')}
                        </p>
                      </div>

                      <div>
                        <Label className="text-slate-400 text-xs">Original Data</Label>
                        <pre className="mt-1 p-3 bg-slate-800 rounded-lg text-xs text-slate-300 overflow-auto max-h-32">
                          {JSON.stringify(selectedItem.originalData, null, 2)}
                        </pre>
                      </div>

                      {selectedItem.suggestedFix && (
                        <div>
                          <Label className="text-slate-400 text-xs flex items-center gap-1">
                            <Sparkles className="h-3 w-3 text-purple-400" />
                            AI Suggested Fix
                          </Label>
                          <pre className="mt-1 p-3 bg-emerald-900/20 border border-emerald-500/30 rounded-lg text-xs text-emerald-300 overflow-auto max-h-32">
                            {JSON.stringify(selectedItem.suggestedFix, null, 2)}
                          </pre>
                        </div>
                      )}

                      <div className="flex gap-2 pt-4 border-t border-slate-700">
                        {selectedItem.status === 'pending' && (
                          <>
                            {!selectedItem.suggestedFix && (
                              <Button
                                size="sm"
                                variant="outline"
                                className="flex-1 border-purple-500 text-purple-400"
                                onClick={() => handleGenerateAISuggestion(selectedItem.id)}
                                disabled={isLoading}
                              >
                                <Sparkles className="h-4 w-4 mr-1" />
                                AI Suggest
                              </Button>
                            )}
                            {selectedItem.suggestedFix && (
                              <Button
                                size="sm"
                                className="flex-1 bg-emerald-600 hover:bg-emerald-700"
                                onClick={() => handleApproveItem(selectedItem.id)}
                              >
                                <Check className="h-4 w-4 mr-1" />
                                Approve
                              </Button>
                            )}
                            <Button
                              size="sm"
                              variant="outline"
                              className="flex-1 border-red-500 text-red-400"
                              onClick={() => handleRejectItem(selectedItem.id)}
                            >
                              <X className="h-4 w-4 mr-1" />
                              Reject
                            </Button>
                          </>
                        )}
                      </div>
                    </div>
                  ) : (
                    <div className="text-center py-12 text-slate-500">
                      <Eye className="h-12 w-12 mx-auto mb-4 opacity-50" />
                      <p>Select an item to view details</p>
                    </div>
                  )}
                </CardContent>
              </Card>
            </div>
          </div>
        </TabsContent>

        {/* Batch AI Fixer Tab */}
        <TabsContent value="batch" className="space-y-6">
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
            {/* Job Controls */}
            <Card className="bg-slate-900/50 border-slate-700">
              <CardHeader>
                <CardTitle className="text-white flex items-center gap-2">
                  <Sparkles className="h-5 w-5 text-purple-400" />
                  Batch AI Fix Jobs
                </CardTitle>
                <CardDescription>
                  Run AI-powered batch corrections on your data
                </CardDescription>
              </CardHeader>
              <CardContent className="space-y-4">
                <div className="grid gap-3">
                  <Button
                    variant="outline"
                    className="justify-start h-auto py-4 border-slate-700 hover:bg-slate-800"
                    onClick={() => handleRunBatchAIFix('smiles_validation')}
                    disabled={isLoading}
                  >
                    <div className="flex items-center gap-3">
                      <div className="p-2 bg-emerald-500/20 rounded-lg">
                        <Beaker className="h-5 w-5 text-emerald-400" />
                      </div>
                      <div className="text-left">
                        <p className="font-medium text-white">SMILES Validation</p>
                        <p className="text-sm text-slate-400">
                          Validate and fix invalid SMILES strings
                        </p>
                      </div>
                    </div>
                  </Button>

                  <Button
                    variant="outline"
                    className="justify-start h-auto py-4 border-slate-700 hover:bg-slate-800"
                    onClick={() => handleRunBatchAIFix('target_normalization')}
                    disabled={isLoading}
                  >
                    <div className="flex items-center gap-3">
                      <div className="p-2 bg-blue-500/20 rounded-lg">
                        <FlaskConical className="h-5 w-5 text-blue-400" />
                      </div>
                      <div className="text-left">
                        <p className="font-medium text-white">Target Normalization</p>
                        <p className="text-sm text-slate-400">
                          Normalize target names to canonical forms
                        </p>
                      </div>
                    </div>
                  </Button>

                  <Button
                    variant="outline"
                    className="justify-start h-auto py-4 border-slate-700 hover:bg-slate-800"
                    onClick={() => handleRunBatchAIFix('embedding_generation')}
                    disabled={isLoading}
                  >
                    <div className="flex items-center gap-3">
                      <div className="p-2 bg-purple-500/20 rounded-lg">
                        <Database className="h-5 w-5 text-purple-400" />
                      </div>
                      <div className="text-left">
                        <p className="font-medium text-white">Embedding Generation</p>
                        <p className="text-sm text-slate-400">
                          Generate missing vector embeddings
                        </p>
                      </div>
                    </div>
                  </Button>
                </div>
              </CardContent>
            </Card>

            {/* Active Jobs */}
            <Card className="bg-slate-900/50 border-slate-700">
              <CardHeader>
                <CardTitle className="text-white">Active Jobs</CardTitle>
              </CardHeader>
              <CardContent>
                <ScrollArea className="h-[300px]">
                  <div className="space-y-3">
                    {batchJobs.length === 0 ? (
                      <div className="text-center py-8 text-slate-500">
                        <Play className="h-8 w-8 mx-auto mb-2 opacity-50" />
                        <p>No batch jobs running</p>
                      </div>
                    ) : (
                      batchJobs.map((job) => (
                        <div
                          key={job.id}
                          className="p-4 bg-slate-800/50 rounded-lg border border-slate-700"
                        >
                          <div className="flex items-center justify-between mb-2">
                            <span className="text-white font-medium">
                              {job.type.replace(/_/g, ' ')}
                            </span>
                            <Badge className={cn(
                              job.status === 'running' && "bg-blue-500/20 text-blue-400",
                              job.status === 'completed' && "bg-emerald-500/20 text-emerald-400",
                              job.status === 'failed' && "bg-red-500/20 text-red-400",
                              job.status === 'pending' && "bg-amber-500/20 text-amber-400"
                            )}>
                              {job.status}
                            </Badge>
                          </div>
                          <Progress
                            value={(job.processedItems / job.totalItems) * 100}
                            className="h-2 mb-2"
                          />
                          <div className="flex justify-between text-xs text-slate-400">
                            <span>{job.processedItems} / {job.totalItems} items</span>
                            <span className="text-emerald-400">
                              {job.successCount} success
                            </span>
                          </div>
                        </div>
                      ))
                    )}
                  </div>
                </ScrollArea>
              </CardContent>
            </Card>
          </div>
        </TabsContent>

        {/* Data Health Tab */}
        <TabsContent value="health" className="space-y-6">
          <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
            {dataHealth.map((health) => (
              <Card key={health.sourceTable} className="bg-slate-900/50 border-slate-700">
                <CardHeader>
                  <CardTitle className="text-white flex items-center gap-2">
                    <Database className="h-5 w-5 text-slate-400" />
                    {health.sourceTable.replace(/_/g, ' ')}
                  </CardTitle>
                  <CardDescription>
                    {health.totalCount.toLocaleString()} total records
                  </CardDescription>
                </CardHeader>
                <CardContent className="space-y-4">
                  <div>
                    <div className="flex justify-between text-sm mb-1">
                      <span className="text-slate-400">Normalization Rate</span>
                      <span className={getHealthColor(health.normalizationRate)}>
                        {health.normalizationRate.toFixed(1)}%
                      </span>
                    </div>
                    <Progress value={health.normalizationRate} className="h-2" />
                    <p className="text-xs text-slate-500 mt-1">
                      {health.normalizedCount.toLocaleString()} / {health.totalCount.toLocaleString()} normalized
                    </p>
                  </div>

                  <div>
                    <div className="flex justify-between text-sm mb-1">
                      <span className="text-slate-400">Embedding Rate</span>
                      <span className={getHealthColor(health.embeddingRate)}>
                        {health.embeddingRate.toFixed(1)}%
                      </span>
                    </div>
                    <Progress value={health.embeddingRate} className="h-2" />
                    <p className="text-xs text-slate-500 mt-1">
                      {health.embeddedCount.toLocaleString()} / {health.totalCount.toLocaleString()} embedded
                    </p>
                  </div>
                </CardContent>
              </Card>
            ))}
          </div>
        </TabsContent>
      </Tabs>
    </div>
  );
}

export default RefinementHubPage;
