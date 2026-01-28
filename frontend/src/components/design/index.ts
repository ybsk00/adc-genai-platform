/**
 * Design Components Index
 * ADC 설계 플랫폼 UI 컴포넌트 모음
 */

// Layout
export { TrinityLayout, TrinityTabsLayout } from './TrinityLayout'

// Visualization
export { MoleculeViewer2D, MoleculeDiffViewer, useRDKit, type ScaffoldHighlight, type MoleculeProperties } from './MoleculeViewer2D'
export { SuccessRadarChart, createDefaultMetrics, type RadarMetric } from './SuccessRadarChart'
export { SAScoreGauge, CompactSAScore } from './SAScoreGauge'

// Evidence & References
export { EvidenceHub, type LiteratureReference, type GoldenSetReference } from './EvidenceHub'

// Audit & Compliance
export { DigitalSealBadge, CompactSealBadge } from './DigitalSealBadge'
export {
  AuditorFeedback,
  CompactAuditBadge,
  type AuditCheckItem,
  type RedesignRequest
} from './AuditorFeedback'

// Optimization Metrics
export {
  DeltaMetricDisplay,
  CompactDelta,
  createOptimizationMetrics,
  type MetricDelta
} from './DeltaMetricDisplay'

// Self-Healing
export { HealerAnimation, HealerStatusIndicator } from './HealerAnimation'

// CMC & Sourcing
export { BOMTable, BOMSummary, type ReagentItem } from './BOMTable'

// UniProt & Antigen Search (Phase 3)
export { AntigenAutocomplete, AntigenBadge } from './AntigenAutocomplete'

// Environment Lock (Phase 3)
export { EnvironmentLock, EnvironmentBadge } from './EnvironmentLock'

// Note: MoleculeViewer2D now uses RDKit.js WASM (upgraded from smiles-drawer)

// Report Export (Phase 4)
export { ReportExport } from './ReportExport'

// Settings (Phase 4 - Theme)
export { SettingsPanel, ThemeToggle } from './SettingsPanel'

// Existing Components
export { LiveAgentConsole } from './LiveAgentConsole'
export { DesignProgressTimeline } from './DesignProgressTimeline'
