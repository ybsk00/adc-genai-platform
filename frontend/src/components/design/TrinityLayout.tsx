/**
 * Trinity Layout - 3단 분할 레이아웃
 * Left: Input Area | Center: Visualization | Right: Evidence Hub
 */
import type { ReactNode } from 'react'
import { cn } from '@/lib/utils'

interface TrinityLayoutProps {
  leftPanel: ReactNode
  centerPanel: ReactNode
  rightPanel: ReactNode
  className?: string
}

export function TrinityLayout({
  leftPanel,
  centerPanel,
  rightPanel,
  className
}: TrinityLayoutProps) {
  return (
    <div className={cn('grid grid-cols-12 gap-4 h-full', className)}>
      {/* Left Panel - Input Area (3 cols on large screens) */}
      <div className="col-span-12 lg:col-span-3 space-y-4 overflow-y-auto">
        {leftPanel}
      </div>

      {/* Center Panel - Visualization Area (5 cols on large screens) */}
      <div className="col-span-12 lg:col-span-5 space-y-4 overflow-y-auto">
        {centerPanel}
      </div>

      {/* Right Panel - Evidence Hub (4 cols on large screens) */}
      <div className="col-span-12 lg:col-span-4 space-y-4 overflow-y-auto">
        {rightPanel}
      </div>
    </div>
  )
}

/**
 * Mobile-friendly Trinity Layout with tabs
 */
import { useState } from 'react'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import { Settings, Eye, BookOpen } from 'lucide-react'

interface TrinityTabsLayoutProps {
  leftPanel: ReactNode
  centerPanel: ReactNode
  rightPanel: ReactNode
  leftLabel?: string
  centerLabel?: string
  rightLabel?: string
}

export function TrinityTabsLayout({
  leftPanel,
  centerPanel,
  rightPanel,
  leftLabel = 'Input',
  centerLabel = 'Visualization',
  rightLabel = 'Evidence'
}: TrinityTabsLayoutProps) {
  const [activeTab, setActiveTab] = useState('center')

  return (
    <>
      {/* Desktop: Grid Layout */}
      <div className="hidden lg:grid grid-cols-12 gap-4 h-full">
        <div className="col-span-3 space-y-4 overflow-y-auto">
          {leftPanel}
        </div>
        <div className="col-span-5 space-y-4 overflow-y-auto">
          {centerPanel}
        </div>
        <div className="col-span-4 space-y-4 overflow-y-auto">
          {rightPanel}
        </div>
      </div>

      {/* Mobile: Tabs Layout */}
      <div className="lg:hidden">
        <Tabs value={activeTab} onValueChange={setActiveTab}>
          <TabsList className="grid w-full grid-cols-3">
            <TabsTrigger value="left" className="flex items-center gap-1">
              <Settings className="w-4 h-4" />
              <span className="hidden sm:inline">{leftLabel}</span>
            </TabsTrigger>
            <TabsTrigger value="center" className="flex items-center gap-1">
              <Eye className="w-4 h-4" />
              <span className="hidden sm:inline">{centerLabel}</span>
            </TabsTrigger>
            <TabsTrigger value="right" className="flex items-center gap-1">
              <BookOpen className="w-4 h-4" />
              <span className="hidden sm:inline">{rightLabel}</span>
            </TabsTrigger>
          </TabsList>
          <TabsContent value="left" className="mt-4">
            {leftPanel}
          </TabsContent>
          <TabsContent value="center" className="mt-4">
            {centerPanel}
          </TabsContent>
          <TabsContent value="right" className="mt-4">
            {rightPanel}
          </TabsContent>
        </Tabs>
      </div>
    </>
  )
}
