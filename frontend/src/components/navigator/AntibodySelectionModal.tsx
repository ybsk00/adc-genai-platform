
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
  DialogDescription,
} from '@/components/ui/dialog';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { ScrollArea } from '@/components/ui/scroll-area';
import { CheckCircle2, AlertTriangle } from 'lucide-react';
import { cn } from '@/lib/utils';

interface AntibodyCandidate {
  id: string;
  name: string;
  target_protein: string;
  clinical_score: number;
  match_confidence: number;
  isotype?: string;
  full_spec?: any;
}

interface AntibodySelectionModalProps {
  open: boolean;
  candidates: AntibodyCandidate[];
  onSelect: (antibodyId: string) => void;
  isLoading: boolean;
}

export function AntibodySelectionModal({
  open,
  candidates,
  onSelect,
  isLoading,
}: AntibodySelectionModalProps) {
  return (
    <Dialog open={open}>
      <DialogContent className="sm:max-w-[600px] bg-slate-900 border-slate-700 text-white">
        <DialogHeader>
          <DialogTitle className="text-xl font-bold flex items-center gap-2">
            <CheckCircle2 className="h-6 w-6 text-violet-400" />
            Select Antibody Candidate
          </DialogTitle>
          <DialogDescription className="text-slate-400">
            AI has identified potential antibody candidates. Please select one to proceed with the design.
          </DialogDescription>
        </DialogHeader>

        <ScrollArea className="h-[400px] pr-4 mt-4">
          <div className="space-y-3">
            {candidates.map((ab) => (
              <div
                key={ab.id}
                className={cn(
                  "p-4 rounded-lg border border-slate-700 bg-slate-800/50 hover:bg-slate-800 transition-colors cursor-pointer group relative overflow-hidden",
                  "hover:border-violet-500/50"
                )}
                onClick={() => !isLoading && onSelect(ab.id)}
              >
                <div className="flex justify-between items-start mb-2">
                  <div>
                    <h4 className="font-semibold text-white group-hover:text-violet-300 transition-colors">
                      {ab.name}
                    </h4>
                    <div className="flex gap-2 mt-1">
                      <Badge variant="secondary" className="bg-slate-700 text-slate-300">
                        {ab.target_protein}
                      </Badge>
                      {ab.isotype && (
                        <Badge variant="outline" className="border-slate-600 text-slate-400">
                          {ab.isotype}
                        </Badge>
                      )}
                    </div>
                  </div>
                  <div className="text-right">
                    <div className="text-sm font-medium text-emerald-400">
                      Match: {(ab.match_confidence * 100).toFixed(0)}%
                    </div>
                    <div className="text-xs text-slate-500 mt-1">
                      Clinical: {(ab.clinical_score * 100).toFixed(0)}
                    </div>
                  </div>
                </div>
                
                {/* Warning for gpNMB (Hardcoded visual cue for demo) */}
                {(ab.target_protein.includes("gpNMB") || ab.target_protein.includes("GPNMB")) && (
                  <div className="mt-2 flex items-center gap-2 text-xs text-amber-400 bg-amber-950/30 p-2 rounded">
                    <AlertTriangle className="h-3 w-3" />
                    High risk of antigen shedding (ID 5552)
                  </div>
                )}
              </div>
            ))}
          </div>
        </ScrollArea>

        <div className="flex justify-end pt-4">
          <Button variant="ghost" disabled>
            Cancel
          </Button>
        </div>
      </DialogContent>
    </Dialog>
  );
}
