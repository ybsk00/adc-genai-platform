/**
 * Navigator Result Component
 * One-Click ADC Navigator - Final Result Display
 *
 * Navigator 파이프라인 완료 후 결과 표시
 */

import { motion } from 'framer-motion';
import {
  CheckCircle2,
  Target,
  Link,
  Pill,
  Award,
  Beaker,
  Download,
  Share2,
  FileText,
  AlertTriangle,
} from 'lucide-react';
import { cn } from '@/lib/utils';
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from '@/components/ui/card';
import { Badge } from '@/components/ui/badge';
import { Button } from '@/components/ui/button';
import { Separator } from '@/components/ui/separator';
import { VirtualTrialChart } from './VirtualTrialChart';

interface AntibodyCandidate {
  id: string;
  name: string;
  target_protein: string;
  clinical_score: number;
  match_confidence: number;
}

interface GoldenCombination {
  antibody: string;
  linker: {
    type: string;
    smiles: string;
  };
  payload: {
    class: string;
    smiles: string;
  };
  dar: number;
  historical_orr: number | null;
}

interface VirtualTrial {
  predicted_orr: number;
  predicted_pfs_months: number;
  predicted_os_months: number;
  pk_data: Array<{ time_hours: number; concentration: number; free_payload: number }>;
  tumor_data: Array<{ day: number; treated: number; control: number }>;
  confidence: number;
}

interface NavigatorResultProps {
  sessionId: string;
  diseaseName: string;
  targetProtein: string | null;
  antibodyCandidates: AntibodyCandidate[];
  goldenCombination: GoldenCombination;
  physicsVerified: boolean;
  virtualTrial: VirtualTrial;
  executionTime: number;
  warnings?: string[];  // FIXED
  dataQualityScore?: number;
  onExportReport?: () => void;
  onShare?: () => void;
}

export function NavigatorResult({
  sessionId,
  diseaseName,
  targetProtein,
  antibodyCandidates,
  goldenCombination,
  physicsVerified,
  virtualTrial,
  executionTime,
  warnings = [],
  dataQualityScore,
  onExportReport,
  onShare,
}: NavigatorResultProps) {
  return (
    <motion.div
      initial={{ opacity: 0, y: 20 }}
      animate={{ opacity: 1, y: 0 }}
      className="space-y-6"
    >
      {/* Header */}
      <div className="text-center">
        <motion.div
          initial={{ scale: 0 }}
          animate={{ scale: 1 }}
          transition={{ type: 'spring', delay: 0.2 }}
          className="inline-flex items-center justify-center w-16 h-16 rounded-full bg-gradient-to-r from-emerald-500 to-teal-500 mb-4"
        >
          <CheckCircle2 className="h-8 w-8 text-white" />
        </motion.div>
        <h2 className="text-2xl font-bold text-white">ADC Design Complete!</h2>
        <p className="text-slate-400 mt-1">
          One-Click Navigator completed in{' '}
          <span className="text-violet-400 font-medium">{Number(executionTime ?? 0).toFixed(1)}s</span>
        </p>

        {/* Physics Verified Badge */}
        {physicsVerified && (
          <motion.div
            initial={{ opacity: 0, scale: 0.8 }}
            animate={{ opacity: 1, scale: 1 }}
            transition={{ delay: 0.4 }}
            className="mt-4 inline-flex items-center gap-2 px-4 py-2 bg-gradient-to-r from-violet-600/20 to-indigo-600/20 border border-violet-500/30 rounded-full"
          >
            <Award className="h-5 w-5 text-violet-400" />
            <span className="text-violet-300 font-medium">Physics Verified</span>
          </motion.div>
        )}
      </div>

      {/* FIXED: Warnings Display */}
      {warnings.length > 0 && (
        <motion.div
          initial={{ opacity: 0, y: -10 }}
          animate={{ opacity: 1, y: 0 }}
          className="bg-amber-900/20 border border-amber-500/30 rounded-lg p-4"
        >
          <div className="flex items-start gap-3">
            <AlertTriangle className="h-5 w-5 text-amber-400 flex-shrink-0 mt-0.5" />
            <div>
              <h4 className="text-amber-300 font-medium">Design Warnings</h4>
              <ul className="mt-2 space-y-1">
                {warnings.map((warning, idx) => (
                  <li key={idx} className="text-sm text-amber-200/80">
                    • {warning}
                  </li>
                ))}
              </ul>
            </div>
          </div>
        </motion.div>
      )}

      {/* FIXED: Data Quality Indicator */}
      {dataQualityScore !== undefined && (
        <div className="flex justify-end">
          <Badge 
            variant="outline" 
            className={cn(
              dataQualityScore > 0.7 ? 'text-emerald-400 border-emerald-400/50' :
              dataQualityScore > 0.4 ? 'text-amber-400 border-amber-400/50' :
              'text-red-400 border-red-400/50'
            )}
          >
            Data Quality: {(dataQualityScore * 100).toFixed(0)}%
          </Badge>
        </div>
      )}

      {/* Main Results Grid */}
      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        {/* Target & Antibody */}
        <Card className="bg-slate-900/50 border-slate-700">
          <CardHeader>
            <CardTitle className="text-white flex items-center gap-2">
              <Target className="h-5 w-5 text-violet-400" />
              Target & Antibody
            </CardTitle>
            <CardDescription>Step 1 Results</CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <div>
              <p className="text-sm text-slate-500">Disease</p>
              <p className="text-lg font-medium text-white">{diseaseName}</p>
            </div>
            <div>
              <p className="text-sm text-slate-500">Primary Target</p>
              <Badge className="mt-1 bg-violet-500/20 text-violet-400">
                {targetProtein || 'Unknown'}
              </Badge>
            </div>
            <Separator className="bg-slate-700" />
            <div>
              <p className="text-sm text-slate-500 mb-2">Top Antibody Candidates</p>
              <div className="space-y-2">
                {antibodyCandidates.slice(0, 3).map((ab, index) => (
                  <div
                    key={ab.id}
                    className={cn(
                      'p-2 rounded-lg',
                      index === 0 ? 'bg-emerald-900/30 border border-emerald-700/50' : 'bg-slate-800/50'
                    )}
                  >
                    <div className="flex items-center justify-between">
                      <span className={cn('font-medium', index === 0 ? 'text-emerald-300' : 'text-slate-300')}>
                        {ab.name}
                      </span>
                      <span className="text-xs text-slate-500">
                        {(Number(ab.match_confidence ?? 0) * 100).toFixed(0)}% match
                      </span>
                    </div>
                  </div>
                ))}
              </div>
            </div>
          </CardContent>
        </Card>

        {/* Golden Combination */}
        <Card className="bg-slate-900/50 border-slate-700">
          <CardHeader>
            <CardTitle className="text-white flex items-center gap-2">
              <Beaker className="h-5 w-5 text-amber-400" />
              Golden Combination
            </CardTitle>
            <CardDescription>Step 2 Results</CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <div>
              <p className="text-sm text-slate-500">Selected Antibody</p>
              <p className="text-lg font-medium text-white">{goldenCombination.antibody}</p>
            </div>

            <div className="grid grid-cols-2 gap-4">
              <div>
                <p className="text-sm text-slate-500 flex items-center gap-1">
                  <Link className="h-3 w-3" />
                  Linker
                </p>
                <Badge
                  variant="outline"
                  className={cn(
                    'mt-1',
                    goldenCombination.linker.type === 'cleavable'
                      ? 'text-emerald-400 border-emerald-400/50'
                      : 'text-amber-400 border-amber-400/50'
                  )}
                >
                  {goldenCombination.linker.type}
                </Badge>
              </div>
              <div>
                <p className="text-sm text-slate-500 flex items-center gap-1">
                  <Pill className="h-3 w-3" />
                  Payload
                </p>
                <Badge variant="outline" className="mt-1 text-violet-400 border-violet-400/50">
                  {goldenCombination.payload.class}
                </Badge>
              </div>
            </div>

            <div className="grid grid-cols-2 gap-4">
              <div className="p-3 bg-slate-800/50 rounded-lg text-center">
                <p className="text-sm text-slate-500">DAR</p>
                <p className="text-2xl font-bold text-white">{goldenCombination.dar}</p>
              </div>
              <div className="p-3 bg-slate-800/50 rounded-lg text-center">
                <p className="text-sm text-slate-500">Historical ORR</p>
                <p className="text-2xl font-bold text-amber-400">
                  {goldenCombination.historical_orr ? `${goldenCombination.historical_orr}%` : 'N/A'}
                </p>
              </div>
            </div>
          </CardContent>
        </Card>

        {/* Predictions Summary */}
        <Card className="bg-slate-900/50 border-slate-700">
          <CardHeader>
            <CardTitle className="text-white flex items-center gap-2">
              <Award className="h-5 w-5 text-emerald-400" />
              Predicted Outcomes
            </CardTitle>
            <CardDescription>Virtual Trial Results</CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <div className="grid grid-cols-1 gap-3">
              <div className="p-3 bg-gradient-to-r from-emerald-900/30 to-emerald-800/30 rounded-lg border border-emerald-700/50">
                <p className="text-sm text-emerald-300/80">Predicted ORR</p>
                <p className="text-3xl font-bold text-emerald-400">
                  {Number(virtualTrial?.predicted_orr ?? 0).toFixed(1)}%
                </p>
              </div>
              <div className="grid grid-cols-2 gap-3">
                <div className="p-3 bg-slate-800/50 rounded-lg text-center">
                  <p className="text-xs text-slate-500">PFS</p>
                  <p className="text-xl font-bold text-white">
                    {Number(virtualTrial?.predicted_pfs_months ?? 0).toFixed(1)}
                    <span className="text-xs text-slate-500 ml-1">mo</span>
                  </p>
                </div>
                <div className="p-3 bg-slate-800/50 rounded-lg text-center">
                  <p className="text-xs text-slate-500">OS</p>
                  <p className="text-xl font-bold text-white">
                    {Number(virtualTrial?.predicted_os_months ?? 0).toFixed(1)}
                    <span className="text-xs text-slate-500 ml-1">mo</span>
                  </p>
                </div>
              </div>
            </div>

            <div className="flex items-center justify-between p-2 bg-slate-800/30 rounded-lg">
              <span className="text-sm text-slate-400">Confidence</span>
              <div className="flex items-center gap-2">
                <div className="w-24 h-2 bg-slate-700 rounded-full overflow-hidden">
                  <div
                    className="h-full bg-gradient-to-r from-violet-500 to-indigo-500"
                    style={{ width: `${virtualTrial.confidence * 100}%` }}
                  />
                </div>
                <span className="text-sm font-medium text-white">
                  {(Number(virtualTrial?.confidence ?? 0) * 100).toFixed(0)}%
                </span>
              </div>
            </div>
          </CardContent>
        </Card>
      </div>

      {/* Virtual Trial Charts */}
      <VirtualTrialChart
        pkData={virtualTrial.pk_data}
        tumorData={virtualTrial.tumor_data}
        predictedORR={virtualTrial.predicted_orr}
      />

      {/* Action Buttons */}
      <div className="flex justify-center gap-4">
        <Button
          variant="outline"
          className="border-slate-700 hover:bg-slate-800"
          onClick={onExportReport}
        >
          <FileText className="mr-2 h-4 w-4" />
          Export Report
        </Button>
        <Button
          variant="outline"
          className="border-slate-700 hover:bg-slate-800"
          onClick={onShare}
        >
          <Share2 className="mr-2 h-4 w-4" />
          Share Results
        </Button>
        <Button className="bg-gradient-to-r from-violet-600 to-indigo-600 hover:from-violet-700 hover:to-indigo-700">
          <Download className="mr-2 h-4 w-4" />
          Download Full Design
        </Button>
      </div>

      {/* Session Info */}
      <div className="text-center text-sm text-slate-600">
        Session ID: {sessionId}
      </div>
    </motion.div>
  );
}

export default NavigatorResult;
