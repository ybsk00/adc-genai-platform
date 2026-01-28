import { BrowserRouter, Routes, Route } from 'react-router-dom'
import { QueryClient, QueryClientProvider } from '@tanstack/react-query'
import { LandingPage } from './pages/Landing'
import { DashboardLayout, DashboardHome, ADCBuilder, ResultViewer, GoldenSetLibrary, DenovoDesign, LeadOptimization, PreclinicalAudit, CMCSourcing } from './pages/Dashboard'
import { AdminLayout, AdminOverview, UserOperations, DataOperations, AITuning, UploadPage, DesignRunsPage, TotalInventoryLayout } from './pages/Admin'
import TotalDataInventory from './pages/Admin/TotalDataInventory'
import { StagingAreaPage } from './pages/Admin/StagingAreaPage'
import { GoldenSetLibraryPage } from './pages/Admin/GoldenSetLibraryPage'
import { LoginPage, SignupPage } from './pages/Auth'
import VerifyReport from './pages/Verify'
import { Toaster } from '@/components/ui/sonner'
import { Layout } from '@/components/layout/Layout'

const queryClient = new QueryClient()

function App() {
  return (
    <QueryClientProvider client={queryClient}>
      <BrowserRouter>
        <Routes>
          {/* Public Routes with Layout */}
          <Route element={<Layout />}>
            <Route path="/" element={<LandingPage />} />
          </Route>

          {/* Auth Routes (No Layout or separate layout if needed) */}
          <Route path="/login" element={<LoginPage />} />
          <Route path="/signup" element={<SignupPage />} />

          {/* Report Verification (Public - QR Code Scan) */}
          <Route path="/verify/:hash" element={<VerifyReport />} />

          {/* Dashboard Routes (Protected) */}
          <Route path="/dashboard" element={<DashboardLayout />}>
            <Route index element={<DashboardHome />} />
            <Route path="builder" element={<ADCBuilder />} />
            <Route path="denovo" element={<DenovoDesign />} />
            <Route path="optimization" element={<LeadOptimization />} />
            <Route path="audit" element={<PreclinicalAudit />} />
            <Route path="cmc" element={<CMCSourcing />} />
            <Route path="result/:jobId" element={<ResultViewer />} />
            <Route path="library" element={<GoldenSetLibrary />} />
          </Route>

          {/* Admin Routes (Protected + Role Check) */}
          <Route path="/admin" element={<AdminLayout />}>
            <Route index element={<AdminOverview />} />
            <Route path="inventory" element={<TotalDataInventory />} />
            <Route path="adc-inventory" element={<TotalInventoryLayout />} />
            <Route path="staging" element={<StagingAreaPage />} />
            <Route path="goldenset" element={<GoldenSetLibraryPage />} />
            <Route path="user-operations" element={<UserOperations />} />
            <Route path="data-operations" element={<DataOperations />} />
            <Route path="upload" element={<UploadPage />} />
            <Route path="ai-tuning" element={<AITuning />} />
            <Route path="design-runs" element={<DesignRunsPage />} />
          </Route>
        </Routes>
        <Toaster />
      </BrowserRouter>
    </QueryClientProvider>
  )
}

export default App
