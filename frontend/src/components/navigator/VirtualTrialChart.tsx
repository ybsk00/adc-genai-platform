/**
 * Virtual Trial Chart Component
 * One-Click ADC Navigator - PK and Tumor Growth Simulation Charts
 *
 * 가상 임상 시뮬레이션 결과 차트 컴포넌트
 */

import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
  ReferenceLine,
} from 'recharts';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Activity, TrendingDown, Clock, Pill } from 'lucide-react';

interface PKDataPoint {
  time_hours: number;
  concentration: number;
  free_payload: number;
}

interface TumorDataPoint {
  day: number;
  treated: number;
  control: number;
}

interface VirtualTrialChartProps {
  pkData: PKDataPoint[];
  tumorData: TumorDataPoint[];
  halfLife?: number;
  tgi?: number;
  predictedORR?: number;
}

// Custom Tooltip for PK Chart
const PKTooltip = ({ active, payload, label }: any) => {
  if (active && payload && payload.length) {
    return (
      <div className="bg-slate-900 border border-slate-700 rounded-lg p-3 shadow-xl">
        <p className="text-slate-400 text-sm mb-2">Time: {label} hours</p>
        {payload.map((entry: any, index: number) => (
          <p key={index} className="text-sm" style={{ color: entry.color }}>
            {entry.name}: {Number(entry.value ?? 0).toFixed(2)} ng/mL
          </p>
        ))}
      </div>
    );
  }
  return null;
};

// Custom Tooltip for Tumor Chart
const TumorTooltip = ({ active, payload, label }: any) => {
  if (active && payload && payload.length) {
    const treated = payload.find((p: any) => p.dataKey === 'treated')?.value;
    const control = payload.find((p: any) => p.dataKey === 'control')?.value;
    const reduction = control && treated ? ((1 - treated / control) * 100).toFixed(1) : null;

    return (
      <div className="bg-slate-900 border border-slate-700 rounded-lg p-3 shadow-xl">
        <p className="text-slate-400 text-sm mb-2">Day {label}</p>
        {payload.map((entry: any, index: number) => (
          <p key={index} className="text-sm" style={{ color: entry.color }}>
            {entry.name}: {Number(entry.value ?? 0).toFixed(0)} mm³
          </p>
        ))}
        {reduction && (
          <p className="text-emerald-400 text-sm mt-1 font-medium">
            Reduction: {reduction}%
          </p>
        )}
      </div>
    );
  }
  return null;
};

export function VirtualTrialChart({
  pkData,
  tumorData,
  halfLife,
  tgi,
  predictedORR,
}: VirtualTrialChartProps) {
  // 데이터 없으면 렌더링하지 않음
  if ((!pkData || pkData.length === 0) && (!tumorData || tumorData.length === 0)) {
    return null;
  }

  // Calculate TGI from tumor data if not provided
  const calculatedTGI = tgi ?? (() => {
    if (tumorData.length < 2) return null;
    const lastPoint = tumorData[tumorData.length - 1];
    const firstPoint = tumorData[0];
    const controlGrowth = lastPoint.control - firstPoint.control;
    const treatedGrowth = lastPoint.treated - firstPoint.treated;
    if (controlGrowth === 0) return null;
    return ((1 - treatedGrowth / controlGrowth) * 100);
  })();

  // Calculate half-life from PK data if not provided
  const calculatedHalfLife = halfLife ?? (() => {
    if (pkData.length < 2) return null;
    const halfConc = pkData[0].concentration / 2;
    const halfPoint = pkData.find((p) => p.concentration <= halfConc);
    return halfPoint?.time_hours;
  })();

  return (
    <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
      {/* PK Chart */}
      <Card className="bg-slate-900/50 border-slate-700">
        <CardHeader className="pb-2">
          <div className="flex items-center justify-between">
            <CardTitle className="text-white flex items-center gap-2">
              <Activity className="h-5 w-5 text-violet-400" />
              Plasma Concentration
            </CardTitle>
            {calculatedHalfLife && (
              <Badge variant="outline" className="text-violet-400 border-violet-400/50">
                <Clock className="h-3 w-3 mr-1" />
                t½ = {Number(calculatedHalfLife ?? 0).toFixed(0)}h
              </Badge>
            )}
          </div>
          <p className="text-sm text-slate-500">
            Pharmacokinetic simulation over time
          </p>
        </CardHeader>
        <CardContent>
          <div className="h-64" style={{ minWidth: 0, minHeight: 0 }}>
            <ResponsiveContainer width="100%" height="100%" minWidth={200} minHeight={200}>
              <LineChart data={pkData} margin={{ top: 5, right: 20, bottom: 5, left: 0 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="#334155" />
                <XAxis
                  dataKey="time_hours"
                  stroke="#64748b"
                  tick={{ fill: '#94a3b8', fontSize: 12 }}
                  label={{
                    value: 'Time (hours)',
                    position: 'insideBottom',
                    offset: -5,
                    fill: '#64748b',
                    fontSize: 12,
                  }}
                />
                <YAxis
                  stroke="#64748b"
                  tick={{ fill: '#94a3b8', fontSize: 12 }}
                  label={{
                    value: 'Concentration (ng/mL)',
                    angle: -90,
                    position: 'insideLeft',
                    fill: '#64748b',
                    fontSize: 12,
                  }}
                />
                <Tooltip content={<PKTooltip />} />
                <Legend
                  wrapperStyle={{ paddingTop: '10px' }}
                  formatter={(value) => (
                    <span className="text-slate-400 text-sm">{value}</span>
                  )}
                />
                {calculatedHalfLife && (
                  <ReferenceLine
                    x={calculatedHalfLife}
                    stroke="#8b5cf6"
                    strokeDasharray="5 5"
                    label={{
                      value: 't½',
                      fill: '#8b5cf6',
                      fontSize: 12,
                    }}
                  />
                )}
                <Line
                  type="monotone"
                  dataKey="concentration"
                  stroke="#8b5cf6"
                  strokeWidth={2}
                  dot={false}
                  name="ADC Concentration"
                />
                <Line
                  type="monotone"
                  dataKey="free_payload"
                  stroke="#22c55e"
                  strokeWidth={2}
                  dot={false}
                  name="Free Payload"
                />
              </LineChart>
            </ResponsiveContainer>
          </div>
        </CardContent>
      </Card>

      {/* Tumor Growth Chart */}
      <Card className="bg-slate-900/50 border-slate-700">
        <CardHeader className="pb-2">
          <div className="flex items-center justify-between">
            <CardTitle className="text-white flex items-center gap-2">
              <TrendingDown className="h-5 w-5 text-emerald-400" />
              Tumor Volume Response
            </CardTitle>
            {calculatedTGI !== null && (
              <Badge
                variant="outline"
                className={`${
                  calculatedTGI > 50
                    ? 'text-emerald-400 border-emerald-400/50'
                    : 'text-amber-400 border-amber-400/50'
                }`}
              >
                <Pill className="h-3 w-3 mr-1" />
                TGI = {Number(calculatedTGI ?? 0).toFixed(0)}%
              </Badge>
            )}
          </div>
          <p className="text-sm text-slate-500">
            Tumor growth inhibition simulation
          </p>
        </CardHeader>
        <CardContent>
          <div className="h-64" style={{ minWidth: 0, minHeight: 0 }}>
            <ResponsiveContainer width="100%" height="100%" minWidth={200} minHeight={200}>
              <LineChart data={tumorData} margin={{ top: 5, right: 20, bottom: 5, left: 0 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="#334155" />
                <XAxis
                  dataKey="day"
                  stroke="#64748b"
                  tick={{ fill: '#94a3b8', fontSize: 12 }}
                  label={{
                    value: 'Day',
                    position: 'insideBottom',
                    offset: -5,
                    fill: '#64748b',
                    fontSize: 12,
                  }}
                />
                <YAxis
                  stroke="#64748b"
                  tick={{ fill: '#94a3b8', fontSize: 12 }}
                  label={{
                    value: 'Volume (mm³)',
                    angle: -90,
                    position: 'insideLeft',
                    fill: '#64748b',
                    fontSize: 12,
                  }}
                />
                <Tooltip content={<TumorTooltip />} />
                <Legend
                  wrapperStyle={{ paddingTop: '10px' }}
                  formatter={(value) => (
                    <span className="text-slate-400 text-sm">{value}</span>
                  )}
                />
                <Line
                  type="monotone"
                  dataKey="control"
                  stroke="#64748b"
                  strokeWidth={2}
                  strokeDasharray="5 5"
                  dot={false}
                  name="Control Group"
                />
                <Line
                  type="monotone"
                  dataKey="treated"
                  stroke="#f97316"
                  strokeWidth={2}
                  dot={false}
                  name="Treated Group"
                />
              </LineChart>
            </ResponsiveContainer>
          </div>
        </CardContent>
      </Card>

      {/* Summary Cards */}
      {predictedORR !== undefined && (
        <Card className="lg:col-span-2 bg-gradient-to-r from-violet-900/30 to-indigo-900/30 border-violet-700/50">
          <CardContent className="py-4">
            <div className="flex items-center justify-center gap-8">
              <div className="text-center">
                <p className="text-slate-400 text-sm">Predicted ORR</p>
                <p className="text-3xl font-bold text-violet-400">
                  {Number(predictedORR ?? 0).toFixed(1)}%
                </p>
              </div>
              {calculatedTGI !== null && (
                <div className="text-center border-l border-slate-700 pl-8">
                  <p className="text-slate-400 text-sm">Tumor Growth Inhibition</p>
                  <p className="text-3xl font-bold text-emerald-400">
                    {Number(calculatedTGI ?? 0).toFixed(1)}%
                  </p>
                </div>
              )}
              {calculatedHalfLife && (
                <div className="text-center border-l border-slate-700 pl-8">
                  <p className="text-slate-400 text-sm">Plasma Half-life</p>
                  <p className="text-3xl font-bold text-amber-400">
                    {Number(calculatedHalfLife ?? 0).toFixed(0)}h
                  </p>
                </div>
              )}
            </div>
          </CardContent>
        </Card>
      )}
    </div>
  );
}

export default VirtualTrialChart;
