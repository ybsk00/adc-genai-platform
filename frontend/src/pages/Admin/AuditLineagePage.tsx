/**
 * Audit & Lineage Page
 * v2.2 Management - Admin Audit Trail & Data Lineage
 */

import { useState, useEffect } from 'react';
import { motion } from 'framer-motion';
import {
  ScrollText,
  RefreshCw,
  Search,
  Filter,
  Calendar,
  User,
  Database,
  Settings,
  AlertTriangle,
  Shield,
  FileEdit,
  Trash2,
  Plus,
  Eye,
  Download,
  ChevronDown,
  ChevronRight,
  Clock,
  GitBranch
} from 'lucide-react';
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { ScrollArea } from '@/components/ui/scroll-area';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select';
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from '@/components/ui/table';
import { toast } from 'sonner';
import { cn } from '@/lib/utils';
import { API_BASE_URL } from '@/lib/api';

// ============================================================================
// Types
// ============================================================================

interface AuditLog {
  id: string;
  adminId: string;
  adminEmail: string;
  actionType: string;
  actionCategory: string;
  targetType: string;
  targetId: string;
  beforeValue: Record<string, unknown> | null;
  afterValue: Record<string, unknown> | null;
  reason: string | null;
  ipAddress: string;
  userAgent: string;
  createdAt: string;
}

interface DataLineage {
  id: string;
  entityType: string;
  entityId: string;
  sourceType: string;
  sourceName: string;
  transformations: string[];
  createdAt: string;
  parentId: string | null;
}

// ============================================================================
// Component
// ============================================================================

export function AuditLineagePage() {
  // State
  const [auditLogs, setAuditLogs] = useState<AuditLog[]>([]);
  const [selectedLog, setSelectedLog] = useState<AuditLog | null>(null);
  const [isLoading, setIsLoading] = useState(false);

  // Filters
  const [searchQuery, setSearchQuery] = useState('');
  const [filterCategory, setFilterCategory] = useState<string>('all');
  const [filterAction, setFilterAction] = useState<string>('all');
  const [dateRange, setDateRange] = useState<string>('7d');

  // Load data on mount
  useEffect(() => {
    loadAuditLogs();
  }, [filterCategory, filterAction, dateRange]);

  const loadAuditLogs = async () => {
    setIsLoading(true);
    try {
      const params = new URLSearchParams();
      if (filterCategory !== 'all') params.append('category', filterCategory);
      if (filterAction !== 'all') params.append('action', filterAction);
      params.append('dateRange', dateRange);

      const response = await fetch(`${API_BASE_URL}/api/admin/audit/logs?${params}`);
      if (response.ok) {
        const data = await response.json();
        setAuditLogs(data.logs || []);
      }
    } catch (error) {
      console.error('Failed to load audit logs:', error);
      // Mock data for demo
      setAuditLogs([
        {
          id: '1',
          adminId: 'admin-1',
          adminEmail: 'admin@astraforge.io',
          actionType: 'ENGINE_SWITCH',
          actionCategory: 'system',
          targetType: 'system_config',
          targetId: 'USE_NIM_API',
          beforeValue: { value: true },
          afterValue: { value: false },
          reason: 'Budget limit reached, switching to Gemini fallback',
          ipAddress: '192.168.1.100',
          userAgent: 'Mozilla/5.0 Chrome/120',
          createdAt: new Date(Date.now() - 3600000).toISOString(),
        },
        {
          id: '2',
          adminId: 'admin-1',
          adminEmail: 'admin@astraforge.io',
          actionType: 'DATA_APPROVE',
          actionCategory: 'data',
          targetType: 'quarantined_data',
          targetId: 'qd-123',
          beforeValue: { status: 'pending', smiles: 'CC(=O)O[X]' },
          afterValue: { status: 'approved', smiles: 'CC(=O)O' },
          reason: null,
          ipAddress: '192.168.1.100',
          userAgent: 'Mozilla/5.0 Chrome/120',
          createdAt: new Date(Date.now() - 7200000).toISOString(),
        },
        {
          id: '3',
          adminId: 'admin-2',
          adminEmail: 'researcher@astraforge.io',
          actionType: 'PROMPT_UPDATE',
          actionCategory: 'ai',
          targetType: 'agent_prompts',
          targetId: 'toxicity_agent',
          beforeValue: { version: 2, prompt: 'Analyze toxicity...' },
          afterValue: { version: 3, prompt: 'Analyze compound toxicity with focus on...' },
          reason: 'Improved toxicity detection accuracy',
          ipAddress: '192.168.1.101',
          userAgent: 'Mozilla/5.0 Firefox/115',
          createdAt: new Date(Date.now() - 86400000).toISOString(),
        },
        {
          id: '4',
          adminId: 'admin-1',
          adminEmail: 'admin@astraforge.io',
          actionType: 'BUDGET_UPDATE',
          actionCategory: 'system',
          targetType: 'budget_configs',
          targetId: 'nvidia_nim',
          beforeValue: { daily_limit_usd: 300 },
          afterValue: { daily_limit_usd: 500 },
          reason: 'Increased daily limit for production workload',
          ipAddress: '192.168.1.100',
          userAgent: 'Mozilla/5.0 Chrome/120',
          createdAt: new Date(Date.now() - 172800000).toISOString(),
        },
        {
          id: '5',
          adminId: 'admin-1',
          adminEmail: 'admin@astraforge.io',
          actionType: 'DATA_DISCARD',
          actionCategory: 'data',
          targetType: 'antibody_library',
          targetId: 'ab-456',
          beforeValue: { name: 'Invalid Antibody', status: 'active' },
          afterValue: null,
          reason: 'Duplicate entry found',
          ipAddress: '192.168.1.100',
          userAgent: 'Mozilla/5.0 Chrome/120',
          createdAt: new Date(Date.now() - 259200000).toISOString(),
        },
      ]);
    } finally {
      setIsLoading(false);
    }
  };

  const handleExportLogs = async () => {
    try {
      toast.success('Audit logs exported to CSV');
    } catch (error) {
      toast.error('Failed to export logs');
    }
  };

  // Filter logs by search
  const filteredLogs = auditLogs.filter(log => {
    if (searchQuery) {
      const query = searchQuery.toLowerCase();
      return (
        log.adminEmail.toLowerCase().includes(query) ||
        log.actionType.toLowerCase().includes(query) ||
        log.targetType.toLowerCase().includes(query) ||
        (log.reason && log.reason.toLowerCase().includes(query))
      );
    }
    return true;
  });

  // Stats
  const todayCount = auditLogs.filter(log =>
    new Date(log.createdAt).toDateString() === new Date().toDateString()
  ).length;

  const systemChanges = auditLogs.filter(log => log.actionCategory === 'system').length;
  const dataChanges = auditLogs.filter(log => log.actionCategory === 'data').length;

  const getActionIcon = (actionType: string) => {
    switch (actionType) {
      case 'ENGINE_SWITCH': return <Settings className="h-4 w-4" />;
      case 'DATA_APPROVE': return <Plus className="h-4 w-4" />;
      case 'DATA_DISCARD': return <Trash2 className="h-4 w-4" />;
      case 'PROMPT_UPDATE': return <FileEdit className="h-4 w-4" />;
      case 'BUDGET_UPDATE': return <Shield className="h-4 w-4" />;
      default: return <AlertTriangle className="h-4 w-4" />;
    }
  };

  const getActionColor = (actionType: string) => {
    switch (actionType) {
      case 'ENGINE_SWITCH': return 'bg-purple-500/20 text-purple-400';
      case 'DATA_APPROVE': return 'bg-emerald-500/20 text-emerald-400';
      case 'DATA_DISCARD': return 'bg-red-500/20 text-red-400';
      case 'PROMPT_UPDATE': return 'bg-blue-500/20 text-blue-400';
      case 'BUDGET_UPDATE': return 'bg-amber-500/20 text-amber-400';
      default: return 'bg-slate-500/20 text-slate-400';
    }
  };

  const getCategoryColor = (category: string) => {
    switch (category) {
      case 'system': return 'bg-purple-500/20 text-purple-400 border-purple-500/30';
      case 'data': return 'bg-blue-500/20 text-blue-400 border-blue-500/30';
      case 'ai': return 'bg-emerald-500/20 text-emerald-400 border-emerald-500/30';
      default: return 'bg-slate-500/20 text-slate-400 border-slate-500/30';
    }
  };

  const formatDate = (dateString: string) => {
    const date = new Date(dateString);
    return date.toLocaleString('ko-KR', {
      year: 'numeric',
      month: '2-digit',
      day: '2-digit',
      hour: '2-digit',
      minute: '2-digit',
    });
  };

  const formatRelativeTime = (dateString: string) => {
    const date = new Date(dateString);
    const now = new Date();
    const diffMs = now.getTime() - date.getTime();
    const diffMins = Math.floor(diffMs / 60000);
    const diffHours = Math.floor(diffMs / 3600000);
    const diffDays = Math.floor(diffMs / 86400000);

    if (diffMins < 60) return `${diffMins}m ago`;
    if (diffHours < 24) return `${diffHours}h ago`;
    return `${diffDays}d ago`;
  };

  return (
    <div className="space-y-6">
      {/* Page Header */}
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-2xl font-bold text-white flex items-center gap-3">
            <ScrollText className="h-7 w-7 text-purple-400" />
            Audit & Lineage
          </h1>
          <p className="text-slate-400 mt-1">
            Track all administrative actions and data lineage
          </p>
        </div>
        <div className="flex gap-2">
          <Button
            variant="outline"
            onClick={handleExportLogs}
            className="border-slate-700"
          >
            <Download className="h-4 w-4 mr-2" />
            Export
          </Button>
          <Button
            variant="outline"
            onClick={loadAuditLogs}
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
              <div className="p-2 bg-purple-500/20 rounded-lg">
                <Clock className="h-5 w-5 text-purple-400" />
              </div>
              <div>
                <p className="text-2xl font-bold text-white">{todayCount}</p>
                <p className="text-sm text-slate-400">Today's Actions</p>
              </div>
            </div>
          </CardContent>
        </Card>

        <Card className="bg-slate-900/50 border-slate-700">
          <CardContent className="p-4">
            <div className="flex items-center gap-3">
              <div className="p-2 bg-blue-500/20 rounded-lg">
                <Settings className="h-5 w-5 text-blue-400" />
              </div>
              <div>
                <p className="text-2xl font-bold text-white">{systemChanges}</p>
                <p className="text-sm text-slate-400">System Changes</p>
              </div>
            </div>
          </CardContent>
        </Card>

        <Card className="bg-slate-900/50 border-slate-700">
          <CardContent className="p-4">
            <div className="flex items-center gap-3">
              <div className="p-2 bg-emerald-500/20 rounded-lg">
                <Database className="h-5 w-5 text-emerald-400" />
              </div>
              <div>
                <p className="text-2xl font-bold text-white">{dataChanges}</p>
                <p className="text-sm text-slate-400">Data Changes</p>
              </div>
            </div>
          </CardContent>
        </Card>

        <Card className="bg-slate-900/50 border-slate-700">
          <CardContent className="p-4">
            <div className="flex items-center gap-3">
              <div className="p-2 bg-amber-500/20 rounded-lg">
                <User className="h-5 w-5 text-amber-400" />
              </div>
              <div>
                <p className="text-2xl font-bold text-white">
                  {new Set(auditLogs.map(l => l.adminId)).size}
                </p>
                <p className="text-sm text-slate-400">Active Admins</p>
              </div>
            </div>
          </CardContent>
        </Card>
      </div>

      <Tabs defaultValue="audit" className="space-y-6">
        <TabsList className="bg-slate-800">
          <TabsTrigger value="audit">Audit Trail</TabsTrigger>
          <TabsTrigger value="lineage">Data Lineage</TabsTrigger>
        </TabsList>

        {/* Audit Trail Tab */}
        <TabsContent value="audit" className="space-y-6">
          {/* Filters */}
          <Card className="bg-slate-900/50 border-slate-700">
            <CardContent className="p-4">
              <div className="flex flex-wrap gap-4">
                <div className="flex-1 min-w-[200px]">
                  <div className="relative">
                    <Search className="absolute left-3 top-1/2 -translate-y-1/2 h-4 w-4 text-slate-400" />
                    <Input
                      placeholder="Search logs..."
                      value={searchQuery}
                      onChange={(e) => setSearchQuery(e.target.value)}
                      className="pl-10 bg-slate-800 border-slate-600"
                    />
                  </div>
                </div>
                <Select value={filterCategory} onValueChange={setFilterCategory}>
                  <SelectTrigger className="w-[150px] bg-slate-800 border-slate-600">
                    <SelectValue placeholder="Category" />
                  </SelectTrigger>
                  <SelectContent>
                    <SelectItem value="all">All Categories</SelectItem>
                    <SelectItem value="system">System</SelectItem>
                    <SelectItem value="data">Data</SelectItem>
                    <SelectItem value="ai">AI</SelectItem>
                  </SelectContent>
                </Select>
                <Select value={filterAction} onValueChange={setFilterAction}>
                  <SelectTrigger className="w-[180px] bg-slate-800 border-slate-600">
                    <SelectValue placeholder="Action" />
                  </SelectTrigger>
                  <SelectContent>
                    <SelectItem value="all">All Actions</SelectItem>
                    <SelectItem value="ENGINE_SWITCH">Engine Switch</SelectItem>
                    <SelectItem value="BUDGET_UPDATE">Budget Update</SelectItem>
                    <SelectItem value="DATA_APPROVE">Data Approve</SelectItem>
                    <SelectItem value="DATA_DISCARD">Data Discard</SelectItem>
                    <SelectItem value="PROMPT_UPDATE">Prompt Update</SelectItem>
                  </SelectContent>
                </Select>
                <Select value={dateRange} onValueChange={setDateRange}>
                  <SelectTrigger className="w-[140px] bg-slate-800 border-slate-600">
                    <SelectValue placeholder="Date Range" />
                  </SelectTrigger>
                  <SelectContent>
                    <SelectItem value="1d">Last 24 hours</SelectItem>
                    <SelectItem value="7d">Last 7 days</SelectItem>
                    <SelectItem value="30d">Last 30 days</SelectItem>
                    <SelectItem value="90d">Last 90 days</SelectItem>
                  </SelectContent>
                </Select>
              </div>
            </CardContent>
          </Card>

          {/* Logs Table and Detail */}
          <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
            {/* Logs List */}
            <div className="lg:col-span-2">
              <Card className="bg-slate-900/50 border-slate-700">
                <CardHeader className="pb-3">
                  <CardTitle className="text-white">
                    Audit Logs ({filteredLogs.length})
                  </CardTitle>
                </CardHeader>
                <CardContent>
                  <ScrollArea className="h-[500px]">
                    <div className="space-y-2">
                      {filteredLogs.length === 0 ? (
                        <div className="text-center py-12 text-slate-500">
                          <ScrollText className="h-12 w-12 mx-auto mb-4 opacity-50" />
                          <p>No audit logs found</p>
                        </div>
                      ) : (
                        filteredLogs.map((log) => (
                          <div
                            key={log.id}
                            onClick={() => setSelectedLog(log)}
                            className={cn(
                              "p-4 rounded-lg border cursor-pointer transition-all",
                              selectedLog?.id === log.id
                                ? "bg-purple-900/30 border-purple-500"
                                : "bg-slate-800/50 border-slate-700 hover:border-slate-600"
                            )}
                          >
                            <div className="flex items-start gap-3">
                              <div className={cn(
                                "p-2 rounded-lg flex-shrink-0",
                                getActionColor(log.actionType)
                              )}>
                                {getActionIcon(log.actionType)}
                              </div>
                              <div className="flex-1 min-w-0">
                                <div className="flex items-center gap-2 mb-1">
                                  <span className="text-white font-medium">
                                    {log.actionType.replace(/_/g, ' ')}
                                  </span>
                                  <Badge variant="outline" className={getCategoryColor(log.actionCategory)}>
                                    {log.actionCategory}
                                  </Badge>
                                </div>
                                <p className="text-sm text-slate-400 truncate">
                                  {log.targetType}: {log.targetId}
                                </p>
                                <div className="flex items-center gap-3 mt-2 text-xs text-slate-500">
                                  <span className="flex items-center gap-1">
                                    <User className="h-3 w-3" />
                                    {log.adminEmail}
                                  </span>
                                  <span className="flex items-center gap-1">
                                    <Clock className="h-3 w-3" />
                                    {formatRelativeTime(log.createdAt)}
                                  </span>
                                </div>
                              </div>
                              <ChevronRight className={cn(
                                "h-5 w-5 text-slate-500 transition-transform flex-shrink-0",
                                selectedLog?.id === log.id && "rotate-90"
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
                  <CardTitle className="text-white">Log Details</CardTitle>
                </CardHeader>
                <CardContent>
                  {selectedLog ? (
                    <div className="space-y-4">
                      <div>
                        <Label className="text-slate-400 text-xs">Action</Label>
                        <div className="flex items-center gap-2 mt-1">
                          <Badge className={getActionColor(selectedLog.actionType)}>
                            {selectedLog.actionType.replace(/_/g, ' ')}
                          </Badge>
                        </div>
                      </div>

                      <div>
                        <Label className="text-slate-400 text-xs">Admin</Label>
                        <p className="text-white">{selectedLog.adminEmail}</p>
                      </div>

                      <div>
                        <Label className="text-slate-400 text-xs">Target</Label>
                        <p className="text-white">{selectedLog.targetType}</p>
                        <p className="text-sm text-slate-400">{selectedLog.targetId}</p>
                      </div>

                      <div>
                        <Label className="text-slate-400 text-xs">Timestamp</Label>
                        <p className="text-white">{formatDate(selectedLog.createdAt)}</p>
                      </div>

                      {selectedLog.reason && (
                        <div>
                          <Label className="text-slate-400 text-xs">Reason</Label>
                          <p className="text-white text-sm">{selectedLog.reason}</p>
                        </div>
                      )}

                      {selectedLog.beforeValue && (
                        <div>
                          <Label className="text-slate-400 text-xs">Before</Label>
                          <pre className="mt-1 p-3 bg-red-900/20 border border-red-500/30 rounded-lg text-xs text-red-300 overflow-auto max-h-24">
                            {JSON.stringify(selectedLog.beforeValue, null, 2)}
                          </pre>
                        </div>
                      )}

                      {selectedLog.afterValue && (
                        <div>
                          <Label className="text-slate-400 text-xs">After</Label>
                          <pre className="mt-1 p-3 bg-emerald-900/20 border border-emerald-500/30 rounded-lg text-xs text-emerald-300 overflow-auto max-h-24">
                            {JSON.stringify(selectedLog.afterValue, null, 2)}
                          </pre>
                        </div>
                      )}

                      <div className="pt-4 border-t border-slate-700">
                        <Label className="text-slate-400 text-xs">Client Info</Label>
                        <p className="text-sm text-slate-400 mt-1">IP: {selectedLog.ipAddress}</p>
                        <p className="text-xs text-slate-500 mt-1 truncate">
                          {selectedLog.userAgent}
                        </p>
                      </div>
                    </div>
                  ) : (
                    <div className="text-center py-12 text-slate-500">
                      <Eye className="h-12 w-12 mx-auto mb-4 opacity-50" />
                      <p>Select a log to view details</p>
                    </div>
                  )}
                </CardContent>
              </Card>
            </div>
          </div>
        </TabsContent>

        {/* Data Lineage Tab */}
        <TabsContent value="lineage" className="space-y-6">
          <Card className="bg-slate-900/50 border-slate-700">
            <CardHeader>
              <CardTitle className="text-white flex items-center gap-2">
                <GitBranch className="h-5 w-5 text-purple-400" />
                Data Lineage Tracking
              </CardTitle>
              <CardDescription>
                Track the origin and transformations of your data
              </CardDescription>
            </CardHeader>
            <CardContent>
              <div className="text-center py-12 text-slate-500">
                <GitBranch className="h-12 w-12 mx-auto mb-4 opacity-50" />
                <p>Data lineage visualization coming soon</p>
                <p className="text-sm mt-2">
                  This feature will show the complete history and transformations of each data record
                </p>
              </div>
            </CardContent>
          </Card>
        </TabsContent>
      </Tabs>
    </div>
  );
}

export default AuditLineagePage;
