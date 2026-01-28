/**
 * Internationalization (i18n) Module
 * English-only translation system
 */

export type Locale = 'en'

// Default locale
const DEFAULT_LOCALE: Locale = 'en'

// Translation type definition
export interface Translations {
  // Common
  common: {
    loading: string
    error: string
    success: string
    cancel: string
    confirm: string
    save: string
    delete: string
    edit: string
    view: string
    download: string
    export: string
    import: string
    search: string
    filter: string
    reset: string
    close: string
    back: string
    next: string
    previous: string
    submit: string
    yes: string
    no: string
  }

  // Navigation
  nav: {
    dashboard: string
    denovoDesign: string
    leadOptimization: string
    preclinicalAudit: string
    cmcSourcing: string
    knowledgeBase: string
    library: string
    settings: string
    profile: string
    logout: string
  }

  // Design Pages
  design: {
    title: string
    subtitle: string
    targetAntigen: string
    targetIndication: string
    dar: string
    linkerPreference: string
    designGoal: string
    startDesign: string
    designInProgress: string
    designCompleted: string
    candidates: string
    generatedCandidates: string
    score: string
    rank: string
    smiles: string
    metrics: string
    premiumOnly: string
  }

  // Audit
  audit: {
    title: string
    subtitle: string
    candidateName: string
    candidateSmiles: string
    auditChecks: string
    startAudit: string
    auditInProgress: string
    auditCompleted: string
    pass: string
    fail: string
    warning: string
    expertReviewRequired: string
  }

  // CMC
  cmc: {
    title: string
    subtitle: string
    productionScale: string
    timeline: string
    budgetRange: string
    specialRequirements: string
    startAnalysis: string
    synthesisRoutes: string
    suppliers: string
    leadTime: string
    estimatedCost: string
    complexity: string
  }

  // Evidence & References
  evidence: {
    title: string
    literatureReferences: string
    goldenSetReferences: string
    noReferences: string
    viewOnPubMed: string
    similarity: string
    clinicalStatus: string
  }

  // Compliance
  compliance: {
    digitalSeal: string
    auditTrail: string
    exportBundle: string
    verified: string
    notVerified: string
    environmentLock: string
    pythonVersion: string
    libraryVersions: string
    cfrCompliant: string
  }

  // Agent
  agent: {
    console: string
    currentAgent: string
    connected: string
    disconnected: string
    healerActive: string
    redesignLoop: string
    iteration: string
  }

  // Reports
  report: {
    title: string
    exportReport: string
    formatPdf: string
    formatHtml: string
    includeCode: string
    preview: string
    generating: string
  }

  // Settings
  settings: {
    theme: string
    lightMode: string
    darkMode: string
    systemDefault: string
  }
}

// English translations
const en: Translations = {
  common: {
    loading: 'Loading...',
    error: 'Error',
    success: 'Success',
    cancel: 'Cancel',
    confirm: 'Confirm',
    save: 'Save',
    delete: 'Delete',
    edit: 'Edit',
    view: 'View',
    download: 'Download',
    export: 'Export',
    import: 'Import',
    search: 'Search',
    filter: 'Filter',
    reset: 'Reset',
    close: 'Close',
    back: 'Back',
    next: 'Next',
    previous: 'Previous',
    submit: 'Submit',
    yes: 'Yes',
    no: 'No'
  },
  nav: {
    dashboard: 'Dashboard',
    denovoDesign: 'De novo Design',
    leadOptimization: 'Lead Optimization',
    preclinicalAudit: 'Pre-clinical Audit',
    cmcSourcing: 'CMC & Sourcing',
    knowledgeBase: 'Knowledge Base',
    library: 'Library',
    settings: 'Settings',
    profile: 'Profile',
    logout: 'Logout'
  },
  design: {
    title: 'De novo Design',
    subtitle: 'AI-powered ADC candidate generation with multi-agent orchestration',
    targetAntigen: 'Target Antigen',
    targetIndication: 'Target Indication',
    dar: 'Drug-to-Antibody Ratio (DAR)',
    linkerPreference: 'Linker Preference',
    designGoal: 'Design Goal',
    startDesign: 'Start De novo Design',
    designInProgress: 'Design in Progress...',
    designCompleted: 'Design completed successfully!',
    candidates: 'Candidates',
    generatedCandidates: 'Generated Candidates',
    score: 'Score',
    rank: 'Rank',
    smiles: 'SMILES',
    metrics: 'Metrics',
    premiumOnly: 'Premium Only'
  },
  audit: {
    title: 'Pre-clinical Audit',
    subtitle: 'Comprehensive safety and efficacy assessment for ADC candidates',
    candidateName: 'Candidate Name',
    candidateSmiles: 'Candidate SMILES',
    auditChecks: 'Audit Checks',
    startAudit: 'Start Audit',
    auditInProgress: 'Audit in Progress...',
    auditCompleted: 'Audit completed',
    pass: 'Pass',
    fail: 'Fail',
    warning: 'Warning',
    expertReviewRequired: 'Expert review required'
  },
  cmc: {
    title: 'CMC & Sourcing',
    subtitle: 'Manufacturing feasibility and supplier identification',
    productionScale: 'Production Scale',
    timeline: 'Timeline',
    budgetRange: 'Budget Range',
    specialRequirements: 'Special Requirements',
    startAnalysis: 'Start CMC Analysis',
    synthesisRoutes: 'Synthesis Routes',
    suppliers: 'Recommended Suppliers',
    leadTime: 'Lead Time',
    estimatedCost: 'Estimated Cost',
    complexity: 'Complexity'
  },
  evidence: {
    title: 'Evidence Hub',
    literatureReferences: 'Literature References',
    goldenSetReferences: 'Golden Set References',
    noReferences: 'No references available',
    viewOnPubMed: 'View on PubMed',
    similarity: 'Similarity',
    clinicalStatus: 'Clinical Status'
  },
  compliance: {
    digitalSeal: 'Digital Seal',
    auditTrail: 'Audit Trail',
    exportBundle: 'Export Bundle',
    verified: 'Verified',
    notVerified: 'Not Verified',
    environmentLock: 'Environment Lock',
    pythonVersion: 'Python Version',
    libraryVersions: 'Library Versions',
    cfrCompliant: '21 CFR Part 11 Compliant'
  },
  agent: {
    console: 'Agent Console',
    currentAgent: 'Current Agent',
    connected: 'Connected',
    disconnected: 'Disconnected',
    healerActive: 'Self-Healing Active',
    redesignLoop: 'Redesign Loop',
    iteration: 'Iteration'
  },
  report: {
    title: 'Report',
    exportReport: 'Export Report',
    formatPdf: 'PDF Format',
    formatHtml: 'HTML Format',
    includeCode: 'Include Code',
    preview: 'Preview',
    generating: 'Generating report...'
  },
  settings: {
    theme: 'Theme',
    lightMode: 'Light Mode',
    darkMode: 'Dark Mode',
    systemDefault: 'System Default'
  }
}

// Translation data
const translations: Record<Locale, Translations> = { en }

// Get current locale (always English)
export function getLocale(): Locale {
  return DEFAULT_LOCALE
}

// Set locale (no-op for English-only)
export function setLocale(_locale: Locale): void {
  // No-op: English only
}

// Get translations
export function getTranslations(_locale?: Locale): Translations {
  return en
}

// Single key translation
export function t(key: string, _locale?: Locale): string {
  const trans = en
  const keys = key.split('.')

  let value: unknown = trans
  for (const k of keys) {
    if (value && typeof value === 'object' && k in value) {
      value = (value as Record<string, unknown>)[k]
    } else {
      return key // fallback to key
    }
  }

  return typeof value === 'string' ? value : key
}

export { translations }
