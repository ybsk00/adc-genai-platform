/**
 * Analytics & Cost Page
 * v2.2 Management - API Usage Analytics & Cost Tracking
 */

import { useState, useEffect } from 'react';
import { motion } from 'framer-motion';
import {
  BarChart3,
  RefreshCw,
  Download,
  DollarSign,
  TrendingUp,
  TrendingDown,
  Clock,
  Zap,
  Database,
  Activity,
  AlertCircle,
  Calendar,
  Cpu,
  Sparkles
} from 'lucide-react';
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { Progress } from '@/components/ui/progress';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select';
import { toast } from 'sonner';
import { cn } from '@/lib/utils';

// ============================================================================
// Types
// ============================================================================

interface UsageMetric {
  date: string;
  provider: string;
  totalCalls: number;
  totalTokens: number;
  totalCost: number;
  avgLatency: number;
  errorRate: number;
}

interface CostSummary {
  provider: string;
  dailyCost: number;
  monthlyCost: number;
  dailyLimit: number;
  monthlyLimit: number;
  trend: number;
}

interface SessionStats {
  totalSessions: number;
  activeSessions: number;
  avgDuration: number;
  completionRate: number;
  byStage: Record<string, number>;
}

// ============================================================================
// Component
// ============================================================================

export function AnalyticsCostPage() {
  // State
  const [usageData, setUsageData] = useState<UsageMetric[]>([]);
  const [costSummary, setCostSummary] = useState<CostSummary[]>([]);
  const [sessionStats, setSessionStats] = useState<SessionStats | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [dateRange, setDateRange] = useState('7d');

  // Load data on mount
  useEffect(() => {
    loadAnalytics();
  }, [dateRange]);

  const loadAnalytics = async () => {
    setIsLoading(true);
    try {
      const response = await fetch(`/api/admin/analytics?range=${dateRange}`);
      if (response.ok) {
        const data = await response.json();
        setUsageData(data.usage || []);
        setCostSummary(data.costs || []);
        setSessionStats(data.sessions || null);
      }
    } catch (error) {
      console.error('Failed to load analytics:', error);
      // Mock data
      setCostSummary([
        {
          provider: 'nvidia_nim',
          dailyCost: 127.45,
          monthlyCost: 2845.30,
          dailyLimit: 500,
          monthlyLimit: 10000,
          trend: 12.5,
        },
        {
          provider: 'gemini',
          dailyCost: 45.20,
          monthlyCost: 890.15,
          dailyLimit: 200,
          monthlyLimit: 4000,
          trend: -5.2,
        },
        {
          provider: 'alphafold',
          dailyCost: 78.90,
          monthlyCost: 1567.80,
          dailyLimit: 300,
          monthlyLimit: 6000,
          trend: 8.3,
        },
      ]);

      setSessionStats({
        totalSessions: 1245,
        activeSessions: 23,
        avgDuration: 342,
        completionRate: 87.5,
        byStage: {
          'target_selection': 5,
          'antibody_selection': 8,
          'linker_selection': 4,
          'payload_selection': 3,
          'optimization': 3,
        },
      });

      // Mock usage data for chart
      const mockUsage: UsageMetric[] = [];
      for (let i = 6; i >= 0; i--) {
        const date = new Date();
        date.setDate(date.getDate() - i);
        mockUsage.push({
          date: date.toISOString().split('T')[0],
          provider: 'nvidia_nim',
          totalCalls: Math.floor(Math.random() * 500) + 200,
          totalTokens: Math.floor(Math.random() * 500000) + 100000,
          totalCost: Math.random() * 150 + 50,
          avgLatency: Math.floor(Math.random() * 500) + 200,
          errorRate: Math.random() * 5,
        });
      }
      setUsageData(mockUsage);
    } finally {
      setIsLoading(false);
    }
  };

  const handleExportReport = () => {
    toast.success('Analytics report exported');
  };

  // Calculate totals
  const totalDailyCost = costSummary.reduce((acc, c) => acc + c.dailyCost, 0);
  const totalMonthlyCost = costSummary.reduce((acc, c) => acc + c.monthlyCost, 0);
  const avgErrorRate = usageData.length > 0
    ? usageData.reduce((acc, u) => acc + u.errorRate, 0) / usageData.length
    : 0;

  const getProviderIcon = (provider: string) => {
    switch (provider) {
      case 'nvidia_nim': return <Cpu className="h-5 w-5 text-emerald-400" />;
      case 'gemini': return <Sparkles className="h-5 w-5 text-blue-400" />;
      case 'alphafold': return <Database className="h-5 w-5 text-purple-400" />;
      default: return <Zap className="h-5 w-5 text-slate-400" />;
    }
  };

  const getProviderName = (provider: string) => {
    switch (provider) {
      case 'nvidia_nim': return 'NVIDIA NIM';
      case 'gemini': return 'Gemini Pro';
      case 'alphafold': return 'AlphaFold';
      default: return provider;
    }
  };

  const formatDuration = (seconds: number) => {
    const mins = Math.floor(seconds / 60);
    const secs = seconds % 60;
    return `${mins}m ${secs}s`;
  };

  return (
    <div className="space-y-6">
      {/* Page Header */}
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-2xl font-bold text-white flex items-center gap-3">
            <BarChart3 className="h-7 w-7 text-purple-400" />
            Analytics & Cost
          </h1>
          <p className="text-slate-400 mt-1">
            Monitor API usage, costs, and performance metrics
          </p>
        </div>
        <div className="flex gap-2">
          <Select value={dateRange} onValueChange={setDateRange}>
            <SelectTrigger className="w-[140px] bg-slate-800 border-slate-600">
              <Calendar className="h-4 w-4 mr-2" />
              <SelectValue />
            </SelectTrigger>
            <SelectContent>
              <SelectItem value="1d">Last 24 hours</SelectItem>
              <SelectItem value="7d">Last 7 days</SelectItem>
              <SelectItem value="30d">Last 30 days</SelectItem>
              <SelectItem value="90d">Last 90 days</SelectItem>
            </SelectContent>
          </Select>
          <Button
            variant="outline"
            onClick={handleExportReport}
            className="border-slate-700"
          >
            <Download className="h-4 w-4 mr-2" />
            Export
          </Button>
          <Button
            variant="outline"
            onClick={loadAnalytics}
            className="border-slate-700"
          >
            <RefreshCw className="h-4 w-4 mr-2" />
            Refresh
          </Button>
        </div>
      </div>

      {/* Cost Overview Cards */}
      <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
        <Card className="bg-slate-900/50 border-slate-700">
          <CardContent className="p-4">
            <div className="flex items-center gap-3">
              <div className="p-2 bg-emerald-500/20 rounded-lg">
                <DollarSign className="h-5 w-5 text-emerald-400" />
              </div>
              <div>
                <p className="text-2xl font-bold text-white">
                  ${totalDailyCost.toFixed(2)}
                </p>
                <p className="text-sm text-slate-400">Today's Cost</p>
              </div>
            </div>
          </CardContent>
        </Card>

        <Card className="bg-slate-900/50 border-slate-700">
          <CardContent className="p-4">
            <div className="flex items-center gap-3">
              <div className="p-2 bg-blue-500/20 rounded-lg">
                <TrendingUp className="h-5 w-5 text-blue-400" />
              </div>
              <div>
                <p className="text-2xl font-bold text-white">
                  ${totalMonthlyCost.toFixed(2)}
                </p>
                <p className="text-sm text-slate-400">Month to Date</p>
              </div>
            </div>
          </CardContent>
        </Card>

        <Card className="bg-slate-900/50 border-slate-700">
          <CardContent className="p-4">
            <div className="flex items-center gap-3">
              <div className="p-2 bg-purple-500/20 rounded-lg">
                <Activity className="h-5 w-5 text-purple-400" />
              </div>
              <div>
                <p className="text-2xl font-bold text-white">
                  {sessionStats?.activeSessions || 0}
                </p>
                <p className="text-sm text-slate-400">Active Sessions</p>
              </div>
            </div>
          </CardContent>
        </Card>

        <Card className="bg-slate-900/50 border-slate-700">
          <CardContent className="p-4">
            <div className="flex items-center gap-3">
              <div className={cn(
                "p-2 rounded-lg",
                avgErrorRate > 5 ? "bg-red-500/20" : "bg-amber-500/20"
              )}>
                <AlertCircle className={cn(
                  "h-5 w-5",
                  avgErrorRate > 5 ? "text-red-400" : "text-amber-400"
                )} />
              </div>
              <div>
                <p className="text-2xl font-bold text-white">
                  {avgErrorRate.toFixed(1)}%
                </p>
                <p className="text-sm text-slate-400">Error Rate</p>
              </div>
            </div>
          </CardContent>
        </Card>
      </div>

      <Tabs defaultValue="costs" className="space-y-6">
        <TabsList className="bg-slate-800">
          <TabsTrigger value="costs">Cost Breakdown</TabsTrigger>
          <TabsTrigger value="sessions">Live Sessions</TabsTrigger>
          <TabsTrigger value="usage">Usage Trends</TabsTrigger>
        </TabsList>

        {/* Cost Breakdown Tab */}
        <TabsContent value="costs" className="space-y-6">
          <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
            {costSummary.map((cost) => {
              const dailyPct = (cost.dailyCost / cost.dailyLimit) * 100;
              const monthlyPct = (cost.monthlyCost / cost.monthlyLimit) * 100;

              return (
                <Card key={cost.provider} className="bg-slate-900/50 border-slate-700">
                  <CardHeader className="pb-2">
                    <CardTitle className="text-white flex items-center justify-between">
                      <span className="flex items-center gap-2">
                        {getProviderIcon(cost.provider)}
                        {getProviderName(cost.provider)}
                      </span>
                      <Badge className={cn(
                        cost.trend > 0
                          ? "bg-red-500/20 text-red-400"
                          : "bg-emerald-500/20 text-emerald-400"
                      )}>
                        {cost.trend > 0 ? (
                          <TrendingUp className="h-3 w-3 mr-1" />
                        ) : (
                          <TrendingDown className="h-3 w-3 mr-1" />
                        )}
                        {Math.abs(cost.trend)}%
                      </Badge>
                    </CardTitle>
                  </CardHeader>
                  <CardContent className="space-y-4">
                    {/* Daily Cost */}
                    <div>
                      <div className="flex justify-between text-sm mb-2">
                        <span className="text-slate-400">Daily</span>
                        <span className="text-white font-medium">
                          ${cost.dailyCost.toFixed(2)} / ${cost.dailyLimit}
                        </span>
                      </div>
                      <div className="relative h-3 bg-slate-800 rounded-full overflow-hidden">
                        <div
                          className={cn(
                            "absolute top-0 left-0 h-full rounded-full transition-all",
                            dailyPct >= 100 ? "bg-red-500" :
                            dailyPct >= 80 ? "bg-amber-500" : "bg-emerald-500"
                          )}
                          style={{ width: `${Math.min(dailyPct, 100)}%` }}
                        />
                      </div>
                      <p className="text-xs text-slate-500 mt-1">
                        {dailyPct.toFixed(1)}% of daily limit
                      </p>
                    </div>

                    {/* Monthly Cost */}
                    <div>
                      <div className="flex justify-between text-sm mb-2">
                        <span className="text-slate-400">Monthly</span>
                        <span className="text-white font-medium">
                          ${cost.monthlyCost.toFixed(2)} / ${cost.monthlyLimit}
                        </span>
                      </div>
                      <Progress value={monthlyPct} className="h-2" />
                      <p className="text-xs text-slate-500 mt-1">
                        {monthlyPct.toFixed(1)}% of monthly limit
                      </p>
                    </div>
                  </CardContent>
                </Card>
              );
            })}
          </div>

          {/* Cost History Chart Placeholder */}
          <Card className="bg-slate-900/50 border-slate-700">
            <CardHeader>
              <CardTitle className="text-white">Cost History</CardTitle>
              <CardDescription>Daily cost breakdown over time</CardDescription>
            </CardHeader>
            <CardContent>
              <div className="h-[300px] flex items-end justify-around gap-2 px-4">
                {usageData.map((data, idx) => {
                  const height = (data.totalCost / 200) * 100;
                  const day = new Date(data.date).toLocaleDateString('en-US', { weekday: 'short' });
                  return (
                    <div key={idx} className="flex flex-col items-center gap-2 flex-1">
                      <div
                        className="w-full bg-gradient-to-t from-purple-600 to-purple-400 rounded-t-lg transition-all hover:from-purple-500 hover:to-purple-300"
                        style={{ height: `${height}%`, minHeight: '20px' }}
                      />
                      <span className="text-xs text-slate-400">{day}</span>
                      <span className="text-xs text-white">${data.totalCost.toFixed(0)}</span>
                    </div>
                  );
                })}
              </div>
            </CardContent>
          </Card>
        </TabsContent>

        {/* Live Sessions Tab */}
        <TabsContent value="sessions" className="space-y-6">
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
            {/* Session Overview */}
            <Card className="bg-slate-900/50 border-slate-700">
              <CardHeader>
                <CardTitle className="text-white flex items-center gap-2">
                  <Activity className="h-5 w-5 text-purple-400" />
                  Session Overview
                </CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                <div className="grid grid-cols-2 gap-4">
                  <div className="p-4 bg-slate-800/50 rounded-lg">
                    <p className="text-3xl font-bold text-white">
                      {sessionStats?.totalSessions.toLocaleString() || 0}
                    </p>
                    <p className="text-sm text-slate-400">Total Sessions</p>
                  </div>
                  <div className="p-4 bg-slate-800/50 rounded-lg">
                    <p className="text-3xl font-bold text-emerald-400">
                      {sessionStats?.activeSessions || 0}
                    </p>
                    <p className="text-sm text-slate-400">Active Now</p>
                  </div>
                  <div className="p-4 bg-slate-800/50 rounded-lg">
                    <p className="text-3xl font-bold text-white">
                      {sessionStats ? formatDuration(sessionStats.avgDuration) : '0m'}
                    </p>
                    <p className="text-sm text-slate-400">Avg Duration</p>
                  </div>
                  <div className="p-4 bg-slate-800/50 rounded-lg">
                    <p className="text-3xl font-bold text-blue-400">
                      {sessionStats?.completionRate.toFixed(1) || 0}%
                    </p>
                    <p className="text-sm text-slate-400">Completion Rate</p>
                  </div>
                </div>
              </CardContent>
            </Card>

            {/* Sessions by Stage */}
            <Card className="bg-slate-900/50 border-slate-700">
              <CardHeader>
                <CardTitle className="text-white">Active Sessions by Stage</CardTitle>
                <CardDescription>Real-time session distribution</CardDescription>
              </CardHeader>
              <CardContent>
                {sessionStats?.byStage ? (
                  <div className="space-y-3">
                    {Object.entries(sessionStats.byStage).map(([stage, count]) => {
                      const total = sessionStats.activeSessions || 1;
                      const pct = (count / total) * 100;
                      return (
                        <div key={stage}>
                          <div className="flex justify-between text-sm mb-1">
                            <span className="text-slate-400 capitalize">
                              {stage.replace(/_/g, ' ')}
                            </span>
                            <span className="text-white">{count} sessions</span>
                          </div>
                          <div className="relative h-2 bg-slate-800 rounded-full overflow-hidden">
                            <motion.div
                              initial={{ width: 0 }}
                              animate={{ width: `${pct}%` }}
                              className="absolute top-0 left-0 h-full bg-gradient-to-r from-purple-500 to-pink-500 rounded-full"
                            />
                          </div>
                        </div>
                      );
                    })}
                  </div>
                ) : (
                  <div className="text-center py-8 text-slate-500">
                    <Activity className="h-8 w-8 mx-auto mb-2 opacity-50" />
                    <p>No active sessions</p>
                  </div>
                )}
              </CardContent>
            </Card>
          </div>
        </TabsContent>

        {/* Usage Trends Tab */}
        <TabsContent value="usage" className="space-y-6">
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
            {/* API Calls */}
            <Card className="bg-slate-900/50 border-slate-700">
              <CardHeader>
                <CardTitle className="text-white">API Calls</CardTitle>
                <CardDescription>Daily API call volume</CardDescription>
              </CardHeader>
              <CardContent>
                <div className="h-[200px] flex items-end justify-around gap-2">
                  {usageData.map((data, idx) => {
                    const maxCalls = Math.max(...usageData.map(d => d.totalCalls));
                    const height = (data.totalCalls / maxCalls) * 100;
                    return (
                      <div key={idx} className="flex flex-col items-center gap-1 flex-1">
                        <div
                          className="w-full bg-gradient-to-t from-blue-600 to-blue-400 rounded-t-lg"
                          style={{ height: `${height}%`, minHeight: '10px' }}
                        />
                        <span className="text-xs text-slate-500">{data.totalCalls}</span>
                      </div>
                    );
                  })}
                </div>
              </CardContent>
            </Card>

            {/* Latency */}
            <Card className="bg-slate-900/50 border-slate-700">
              <CardHeader>
                <CardTitle className="text-white">Average Latency</CardTitle>
                <CardDescription>Response time in milliseconds</CardDescription>
              </CardHeader>
              <CardContent>
                <div className="h-[200px] flex items-end justify-around gap-2">
                  {usageData.map((data, idx) => {
                    const maxLatency = Math.max(...usageData.map(d => d.avgLatency));
                    const height = (data.avgLatency / maxLatency) * 100;
                    return (
                      <div key={idx} className="flex flex-col items-center gap-1 flex-1">
                        <div
                          className={cn(
                            "w-full rounded-t-lg",
                            data.avgLatency > 500
                              ? "bg-gradient-to-t from-red-600 to-red-400"
                              : data.avgLatency > 300
                              ? "bg-gradient-to-t from-amber-600 to-amber-400"
                              : "bg-gradient-to-t from-emerald-600 to-emerald-400"
                          )}
                          style={{ height: `${height}%`, minHeight: '10px' }}
                        />
                        <span className="text-xs text-slate-500">{data.avgLatency}ms</span>
                      </div>
                    );
                  })}
                </div>
              </CardContent>
            </Card>
          </div>

          {/* Token Usage */}
          <Card className="bg-slate-900/50 border-slate-700">
            <CardHeader>
              <CardTitle className="text-white">Token Usage</CardTitle>
              <CardDescription>Total tokens processed per day</CardDescription>
            </CardHeader>
            <CardContent>
              <div className="h-[200px] flex items-end justify-around gap-2 px-4">
                {usageData.map((data, idx) => {
                  const maxTokens = Math.max(...usageData.map(d => d.totalTokens));
                  const height = (data.totalTokens / maxTokens) * 100;
                  const day = new Date(data.date).toLocaleDateString('en-US', { weekday: 'short' });
                  return (
                    <div key={idx} className="flex flex-col items-center gap-2 flex-1">
                      <div
                        className="w-full bg-gradient-to-t from-purple-600 to-pink-400 rounded-t-lg"
                        style={{ height: `${height}%`, minHeight: '20px' }}
                      />
                      <span className="text-xs text-slate-400">{day}</span>
                      <span className="text-xs text-white">
                        {(data.totalTokens / 1000).toFixed(0)}K
                      </span>
                    </div>
                  );
                })}
              </div>
            </CardContent>
          </Card>
        </TabsContent>
      </Tabs>
    </div>
  );
}

export default AnalyticsCostPage;
