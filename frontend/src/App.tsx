import { BrowserRouter, Routes, Route } from 'react-router-dom'
import { QueryClient, QueryClientProvider } from '@tanstack/react-query'
import { LandingPage } from './pages/Landing'
import { DashboardLayout, DashboardHome, ADCBuilder, ResultViewer, GoldenSetLibrary } from './pages/Dashboard'
import { AdminLayout, AdminOverview, UserOperations, DataOperations, AITuning } from './pages/Admin'
import { LoginPage, SignupPage } from './pages/Auth'
import { Toaster } from '@/components/ui/sonner'

const queryClient = new QueryClient()

function App() {
  return (
    <QueryClientProvider client={queryClient}>
      <BrowserRouter>
        <Routes>
          {/* Public Routes */}
          <Route path="/" element={<LandingPage />} />
          <Route path="/login" element={<LoginPage />} />
          <Route path="/signup" element={<SignupPage />} />

          {/* Dashboard Routes (Protected) */}
          <Route path="/dashboard" element={<DashboardLayout />}>
            <Route index element={<DashboardHome />} />
            <Route path="builder" element={<ADCBuilder />} />
            <Route path="result/:jobId" element={<ResultViewer />} />
            <Route path="library" element={<GoldenSetLibrary />} />
          </Route>

          {/* Admin Routes (Protected + Role Check) */}
          <Route path="/admin" element={<AdminLayout />}>
            <Route index element={<AdminOverview />} />
            <Route path="users" element={<UserOperations />} />
            <Route path="data" element={<DataOperations />} />
            <Route path="ai" element={<AITuning />} />
          </Route>
        </Routes>
        <Toaster />
      </BrowserRouter>
    </QueryClientProvider>
  )
}

export default App
