/**
 * AI Engine Control Page
 * v2.2 Management - Engine Switch & Emergency Brake
 */

import { useState, useEffect, useCallback } from 'react';
import { motion } from 'framer-motion';
import {
  Cpu,
  Zap,
  AlertTriangle,
  Shield,
  RefreshCw,
  Save,
  History,
  DollarSign,
  Gauge,
  Power,
  Bell,
  ExternalLink
} from 'lucide-react';
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { Switch } from '@/components/ui/switch';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Badge } from '@/components/ui/badge';
import { Progress } from '@/components/ui/progress';
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs';
import { toast } from 'sonner';
import { cn } from '@/lib/utils';
import { API_BASE_URL } from '@/lib/api';

// ============================================================================
// Types
// ============================================================================

interface EngineConfig {
  useNimApi: boolean;
  fallbackEnabled: boolean;
  emergencyStop: boolean;
}

interface BudgetConfig {
  apiProvider: string;
  dailyLimitUsd: number;
  currentUsageUsd: number;
  monthlyLimitUsd: number;
  currentMonthlyUsd: number;
  autoFallbackEnabled: boolean;
  alertThresholdPct: number;
  alertEmail: string;
}

interface AgentPrompt {
  id: string;
  agentName: string;
  version: number;
  systemPrompt: string;
  isActive: boolean;
  createdAt: string;
}

// ============================================================================
// Component
// ============================================================================

export function EngineControlPage() {
  // Engine Config State
  const [engineConfig, setEngineConfig] = useState<EngineConfig>({
    useNimApi: true,
    fallbackEnabled: true,
    emergencyStop: false,
  });

  // Budget State
  const [budgets, setBudgets] = useState<BudgetConfig[]>([
    {
      apiProvider: 'nvidia_nim',
      dailyLimitUsd: 500,
      currentUsageUsd: 0,
      monthlyLimitUsd: 10000,
      currentMonthlyUsd: 0,
      autoFallbackEnabled: true,
      alertThresholdPct: 80,
      alertEmail: '',
    },
    {
      apiProvider: 'gemini',
      dailyLimitUsd: 200,
      currentUsageUsd: 0,
      monthlyLimitUsd: 4000,
      currentMonthlyUsd: 0,
      autoFallbackEnabled: true,
      alertThresholdPct: 80,
      alertEmail: '',
    },
  ]);

  const [isLoading, setIsLoading] = useState(false);
  const [isSaving, setIsSaving] = useState(false);

  // Load initial data
  useEffect(() => {
    loadEngineConfig();
    loadBudgetConfig();
  }, []);

  const loadEngineConfig = async () => {
    try {
      const response = await fetch(`${API_BASE_URL}/api/admin/engine/config`);
      if (response.ok) {
        const data = await response.json();
        setEngineConfig(data);
      }
    } catch (error) {
      console.error('Failed to load engine config:', error);
    }
  };

  const loadBudgetConfig = async () => {
    try {
      const response = await fetch(`${API_BASE_URL}/api/admin/budget/status`);
      if (response.ok) {
        const data = await response.json();
        if (data.budgets) {
          setBudgets(data.budgets);
        }
      }
    } catch (error) {
      console.error('Failed to load budget config:', error);
    }
  };

  // Handlers
  const handleEngineSwitch = async (useNim: boolean) => {
    setIsSaving(true);
    try {
      const response = await fetch(`${API_BASE_URL}/api/admin/engine/switch`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ useNimApi: useNim }),
      });

      if (response.ok) {
        setEngineConfig(prev => ({ ...prev, useNimApi: useNim }));
        toast.success(`Engine switched to ${useNim ? 'NVIDIA NIM' : 'Gemini Fallback'}`);
      } else {
        toast.error('Failed to switch engine');
      }
    } catch (error) {
      toast.error('Failed to switch engine');
    } finally {
      setIsSaving(false);
    }
  };

  const handleEmergencyStop = async () => {
    if (!confirm('Are you sure you want to activate Emergency Stop? This will halt all API calls.')) {
      return;
    }

    setIsSaving(true);
    try {
      const response = await fetch(`${API_BASE_URL}/api/admin/budget/emergency-stop`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ activate: true }),
      });

      if (response.ok) {
        setEngineConfig(prev => ({ ...prev, emergencyStop: true }));
        toast.warning('Emergency Stop activated! All API calls are halted.');
      }
    } catch (error) {
      toast.error('Failed to activate Emergency Stop');
    } finally {
      setIsSaving(false);
    }
  };

  const handleDeactivateEmergencyStop = async () => {
    setIsSaving(true);
    try {
      const response = await fetch(`${API_BASE_URL}/api/admin/budget/emergency-stop`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ activate: false }),
      });

      if (response.ok) {
        setEngineConfig(prev => ({ ...prev, emergencyStop: false }));
        toast.success('Emergency Stop deactivated. Normal operations resumed.');
      }
    } catch (error) {
      toast.error('Failed to deactivate Emergency Stop');
    } finally {
      setIsSaving(false);
    }
  };

  const handleSaveBudget = async (provider: string, config: Partial<BudgetConfig>) => {
    setIsSaving(true);
    try {
      const response = await fetch(`${API_BASE_URL}/api/admin/budget/config`, {
        method: 'PUT',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ provider, ...config }),
      });

      if (response.ok) {
        toast.success('Budget configuration saved');
        loadBudgetConfig();
      }
    } catch (error) {
      toast.error('Failed to save budget config');
    } finally {
      setIsSaving(false);
    }
  };

  // Calculate usage percentage
  const getUsagePercentage = (current: number, limit: number) => {
    if (limit === 0) return 0;
    return Math.min((current / limit) * 100, 100);
  };

  const getUsageColor = (pct: number) => {
    if (pct >= 100) return 'bg-red-500';
    if (pct >= 80) return 'bg-amber-500';
    return 'bg-emerald-500';
  };

  return (
    <div className="space-y-6">
      {/* Page Header */}
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-2xl font-bold text-white flex items-center gap-3">
            <Cpu className="h-7 w-7 text-purple-400" />
            AI Engine Control
          </h1>
          <p className="text-slate-400 mt-1">
            Manage AI engine settings, budget limits, and emergency controls
          </p>
        </div>
        <Button
          variant="outline"
          onClick={() => { loadEngineConfig(); loadBudgetConfig(); }}
          className="border-slate-700"
        >
          <RefreshCw className="h-4 w-4 mr-2" />
          Refresh
        </Button>
      </div>

      {/* Emergency Stop Banner */}
      {engineConfig.emergencyStop && (
        <motion.div
          initial={{ opacity: 0, y: -20 }}
          animate={{ opacity: 1, y: 0 }}
          className="p-4 bg-red-900/50 border border-red-500 rounded-lg flex items-center justify-between"
        >
          <div className="flex items-center gap-3">
            <AlertTriangle className="h-6 w-6 text-red-400 animate-pulse" />
            <div>
              <h3 className="text-red-300 font-semibold">Emergency Stop Active</h3>
              <p className="text-red-400 text-sm">All API calls are currently halted</p>
            </div>
          </div>
          <Button
            variant="outline"
            className="border-red-500 text-red-400 hover:bg-red-900/50"
            onClick={handleDeactivateEmergencyStop}
            disabled={isSaving}
          >
            <Power className="h-4 w-4 mr-2" />
            Deactivate
          </Button>
        </motion.div>
      )}

      <Tabs defaultValue="engine" className="space-y-6">
        <TabsList className="bg-slate-800">
          <TabsTrigger value="engine">Engine Switch</TabsTrigger>
          <TabsTrigger value="budget">Budget Control</TabsTrigger>
          <TabsTrigger value="prompts">Agent Prompts</TabsTrigger>
        </TabsList>

        {/* Engine Switch Tab */}
        <TabsContent value="engine" className="space-y-6">
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
            {/* Current Engine Status */}
            <Card className="bg-slate-900/50 border-slate-700">
              <CardHeader>
                <CardTitle className="text-white flex items-center gap-2">
                  <Zap className="h-5 w-5 text-amber-400" />
                  Current Engine Status
                </CardTitle>
                <CardDescription>Active AI inference engine</CardDescription>
              </CardHeader>
              <CardContent className="space-y-6">
                <div className="flex items-center justify-between p-4 bg-slate-800/50 rounded-lg">
                  <div className="flex items-center gap-3">
                    <div className={cn(
                      "w-12 h-12 rounded-lg flex items-center justify-center",
                      engineConfig.useNimApi
                        ? "bg-emerald-500/20 text-emerald-400"
                        : "bg-blue-500/20 text-blue-400"
                    )}>
                      {engineConfig.useNimApi ? (
                        <Cpu className="h-6 w-6" />
                      ) : (
                        <Zap className="h-6 w-6" />
                      )}
                    </div>
                    <div>
                      <p className="font-semibold text-white">
                        {engineConfig.useNimApi ? 'NVIDIA NIM' : 'Gemini Pro'}
                      </p>
                      <p className="text-sm text-slate-400">
                        {engineConfig.useNimApi
                          ? 'High-performance GPU inference'
                          : 'Fallback mode active'}
                      </p>
                    </div>
                  </div>
                  <Badge className={cn(
                    engineConfig.useNimApi
                      ? "bg-emerald-500/20 text-emerald-400"
                      : "bg-blue-500/20 text-blue-400"
                  )}>
                    Active
                  </Badge>
                </div>

                <div className="space-y-4">
                  <div className="flex items-center justify-between">
                    <div>
                      <Label className="text-white">Use NVIDIA NIM API</Label>
                      <p className="text-sm text-slate-400">
                        Primary inference engine for optimal performance
                      </p>
                    </div>
                    <Switch
                      checked={engineConfig.useNimApi}
                      onCheckedChange={handleEngineSwitch}
                      disabled={isSaving || engineConfig.emergencyStop}
                    />
                  </div>

                  <div className="flex items-center justify-between">
                    <div>
                      <Label className="text-white">Auto-Fallback Enabled</Label>
                      <p className="text-sm text-slate-400">
                        Automatically switch to Gemini when budget exceeded
                      </p>
                    </div>
                    <Switch
                      checked={engineConfig.fallbackEnabled}
                      onCheckedChange={(checked) =>
                        setEngineConfig(prev => ({ ...prev, fallbackEnabled: checked }))
                      }
                      disabled={isSaving}
                    />
                  </div>
                </div>
              </CardContent>
            </Card>

            {/* Emergency Brake */}
            <Card className="bg-slate-900/50 border-slate-700">
              <CardHeader>
                <CardTitle className="text-white flex items-center gap-2">
                  <AlertTriangle className="h-5 w-5 text-red-400" />
                  Emergency Brake
                </CardTitle>
                <CardDescription>Immediate API shutdown controls</CardDescription>
              </CardHeader>
              <CardContent className="space-y-4">
                <div className="p-4 bg-red-900/20 border border-red-500/30 rounded-lg">
                  <p className="text-sm text-red-300 mb-4">
                    Emergency Stop will immediately halt all API calls to prevent
                    unexpected cost overruns. Use with caution.
                  </p>
                  <Button
                    variant="destructive"
                    className="w-full"
                    onClick={handleEmergencyStop}
                    disabled={isSaving || engineConfig.emergencyStop}
                  >
                    <AlertTriangle className="h-4 w-4 mr-2" />
                    {engineConfig.emergencyStop ? 'Emergency Stop Active' : 'Activate Emergency Stop'}
                  </Button>
                </div>

                <div className="space-y-2">
                  <Label className="text-slate-400 text-sm">Quick Actions</Label>
                  <div className="grid grid-cols-2 gap-2">
                    <Button
                      variant="outline"
                      size="sm"
                      className="border-slate-700"
                      onClick={() => handleEngineSwitch(false)}
                      disabled={!engineConfig.useNimApi || isSaving}
                    >
                      Switch to Gemini
                    </Button>
                    <Button
                      variant="outline"
                      size="sm"
                      className="border-slate-700"
                      onClick={() => handleEngineSwitch(true)}
                      disabled={engineConfig.useNimApi || isSaving}
                    >
                      Switch to NIM
                    </Button>
                  </div>
                </div>
              </CardContent>
            </Card>
          </div>
        </TabsContent>

        {/* Budget Control Tab */}
        <TabsContent value="budget" className="space-y-6">
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
            {budgets.map((budget) => {
              const dailyPct = getUsagePercentage(budget.currentUsageUsd, budget.dailyLimitUsd);
              const monthlyPct = getUsagePercentage(budget.currentMonthlyUsd, budget.monthlyLimitUsd);

              return (
                <Card key={budget.apiProvider} className="bg-slate-900/50 border-slate-700">
                  <CardHeader>
                    <CardTitle className="text-white flex items-center justify-between">
                      <span className="flex items-center gap-2">
                        <DollarSign className="h-5 w-5 text-emerald-400" />
                        {budget.apiProvider === 'nvidia_nim' ? 'NVIDIA NIM' : 'Gemini Pro'}
                      </span>
                      {dailyPct >= budget.alertThresholdPct && (
                        <Badge className="bg-amber-500/20 text-amber-400 animate-pulse">
                          <Bell className="h-3 w-3 mr-1" />
                          Alert
                        </Badge>
                      )}
                    </CardTitle>
                  </CardHeader>
                  <CardContent className="space-y-6">
                    {/* Daily Usage */}
                    <div>
                      <div className="flex justify-between text-sm mb-2">
                        <span className="text-slate-400">Daily Usage</span>
                        <span className="text-white font-medium">
                          ${budget.currentUsageUsd.toFixed(2)} / ${budget.dailyLimitUsd.toFixed(2)}
                        </span>
                      </div>
                      <div className="relative">
                        <Progress
                          value={dailyPct}
                          className="h-3"
                        />
                        <div
                          className={cn(
                            "absolute top-0 left-0 h-3 rounded-full transition-all",
                            getUsageColor(dailyPct)
                          )}
                          style={{ width: `${dailyPct}%` }}
                        />
                      </div>
                      <p className="text-xs text-slate-500 mt-1">
                        {dailyPct.toFixed(1)}% used
                        {dailyPct >= budget.alertThresholdPct && ' - Threshold exceeded!'}
                      </p>
                    </div>

                    {/* Monthly Usage */}
                    <div>
                      <div className="flex justify-between text-sm mb-2">
                        <span className="text-slate-400">Monthly Usage</span>
                        <span className="text-white font-medium">
                          ${budget.currentMonthlyUsd.toFixed(2)} / ${budget.monthlyLimitUsd.toFixed(2)}
                        </span>
                      </div>
                      <Progress value={monthlyPct} className="h-2" />
                    </div>

                    {/* Settings */}
                    <div className="space-y-4 pt-4 border-t border-slate-700">
                      <div className="grid grid-cols-2 gap-4">
                        <div>
                          <Label className="text-slate-400 text-xs">Daily Limit ($)</Label>
                          <Input
                            type="number"
                            value={budget.dailyLimitUsd}
                            onChange={(e) => {
                              const value = parseFloat(e.target.value);
                              setBudgets(prev => prev.map(b =>
                                b.apiProvider === budget.apiProvider
                                  ? { ...b, dailyLimitUsd: value }
                                  : b
                              ));
                            }}
                            className="bg-slate-800 border-slate-600 mt-1"
                          />
                        </div>
                        <div>
                          <Label className="text-slate-400 text-xs">Alert at (%)</Label>
                          <Input
                            type="number"
                            value={budget.alertThresholdPct}
                            onChange={(e) => {
                              const value = parseInt(e.target.value);
                              setBudgets(prev => prev.map(b =>
                                b.apiProvider === budget.apiProvider
                                  ? { ...b, alertThresholdPct: value }
                                  : b
                              ));
                            }}
                            className="bg-slate-800 border-slate-600 mt-1"
                          />
                        </div>
                      </div>

                      <div className="flex items-center justify-between">
                        <div>
                          <Label className="text-white text-sm">Auto-Fallback</Label>
                          <p className="text-xs text-slate-500">
                            Switch to Gemini at 100%
                          </p>
                        </div>
                        <Switch
                          checked={budget.autoFallbackEnabled}
                          onCheckedChange={(checked) => {
                            setBudgets(prev => prev.map(b =>
                              b.apiProvider === budget.apiProvider
                                ? { ...b, autoFallbackEnabled: checked }
                                : b
                            ));
                          }}
                        />
                      </div>

                      <Button
                        size="sm"
                        className="w-full"
                        onClick={() => handleSaveBudget(budget.apiProvider, budget)}
                        disabled={isSaving}
                      >
                        <Save className="h-4 w-4 mr-2" />
                        Save Changes
                      </Button>
                    </div>
                  </CardContent>
                </Card>
              );
            })}
          </div>
        </TabsContent>

        {/* Agent Prompts Tab */}
        <TabsContent value="prompts" className="space-y-6">
          <Card className="bg-slate-900/50 border-slate-700">
            <CardHeader>
              <CardTitle className="text-white">Agent Prompt Management</CardTitle>
              <CardDescription>
                Manage system prompts for each AI agent with version control
              </CardDescription>
            </CardHeader>
            <CardContent>
              <div className="text-center py-12 text-slate-500">
                <History className="h-12 w-12 mx-auto mb-4 opacity-50" />
                <p>Prompt management interface coming soon</p>
                <p className="text-sm mt-2">
                  This feature will allow version-controlled prompt editing
                </p>
              </div>
            </CardContent>
          </Card>
        </TabsContent>
      </Tabs>
    </div>
  );
}

export default EngineControlPage;
