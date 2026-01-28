/**
 * SettingsPanel Component
 * Theme settings panel (Dark/Light mode)
 */
import { useState } from 'react'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Label } from '@/components/ui/label'
import {
  Dialog,
  DialogContent,
  DialogDescription,
  DialogHeader,
  DialogTitle,
  DialogTrigger
} from '@/components/ui/dialog'
import {
  Settings,
  Palette,
  Sun,
  Moon,
  Monitor,
  Check
} from 'lucide-react'
import { useTheme } from '@/hooks/useTheme'
import type { Theme } from '@/hooks/useTheme'
import { cn } from '@/lib/utils'

interface SettingsPanelProps {
  compact?: boolean
  className?: string
}

export function SettingsPanel({ compact = false, className }: SettingsPanelProps) {
  const [isOpen, setIsOpen] = useState(false)
  const { theme, setTheme } = useTheme()

  // Theme options
  const themeOptions: { value: Theme; label: string; icon: React.ReactNode }[] = [
    { value: 'light', label: 'Light Mode', icon: <Sun className="w-4 h-4" /> },
    { value: 'dark', label: 'Dark Mode', icon: <Moon className="w-4 h-4" /> },
    { value: 'system', label: 'System Default', icon: <Monitor className="w-4 h-4" /> }
  ]

  // Compact version (icon button)
  if (compact) {
    return (
      <Dialog open={isOpen} onOpenChange={setIsOpen}>
        <DialogTrigger asChild>
          <Button variant="ghost" size="icon" className={className}>
            <Settings className="w-5 h-5" />
          </Button>
        </DialogTrigger>
        <DialogContent className="sm:max-w-md">
          <DialogHeader>
            <DialogTitle className="flex items-center gap-2">
              <Settings className="w-5 h-5" />
              Settings
            </DialogTitle>
            <DialogDescription>
              Theme settings
            </DialogDescription>
          </DialogHeader>

          <SettingsContent
            theme={theme}
            setTheme={setTheme}
            themeOptions={themeOptions}
          />
        </DialogContent>
      </Dialog>
    )
  }

  // Full version (Card)
  return (
    <Card className={className}>
      <CardHeader className="pb-3">
        <CardTitle className="text-sm flex items-center gap-2">
          <Settings className="w-4 h-4" />
          Settings
        </CardTitle>
        <CardDescription className="text-xs">
          Theme settings
        </CardDescription>
      </CardHeader>
      <CardContent>
        <SettingsContent
          theme={theme}
          setTheme={setTheme}
          themeOptions={themeOptions}
        />
      </CardContent>
    </Card>
  )
}

// Settings content component
function SettingsContent({
  theme,
  setTheme,
  themeOptions
}: {
  theme: Theme
  setTheme: (theme: Theme) => void
  themeOptions: { value: Theme; label: string; icon: React.ReactNode }[]
}) {
  return (
    <div className="space-y-4">
      {/* Theme Settings */}
      <div className="space-y-3">
        <Label className="flex items-center gap-2">
          <Palette className="w-4 h-4" />
          Theme
        </Label>
        <div className="grid grid-cols-3 gap-2">
          {themeOptions.map((option) => (
            <button
              key={option.value}
              onClick={() => setTheme(option.value)}
              className={cn(
                'flex flex-col items-center gap-2 p-3 rounded-lg border transition-all',
                theme === option.value
                  ? 'border-blue-500 bg-blue-50 dark:bg-blue-950'
                  : 'border-gray-200 hover:border-gray-300 dark:border-gray-700'
              )}
            >
              {option.icon}
              <span className="text-xs font-medium">{option.label}</span>
              {theme === option.value && (
                <Check className="w-3 h-3 text-blue-500" />
              )}
            </button>
          ))}
        </div>
      </div>
    </div>
  )
}

// Quick Theme Toggle (for header)
export function ThemeToggle({ className }: { className?: string }) {
  const { resolvedTheme, toggleTheme } = useTheme()

  return (
    <Button
      variant="ghost"
      size="icon"
      className={className}
      onClick={toggleTheme}
    >
      {resolvedTheme === 'dark' ? (
        <Sun className="w-4 h-4" />
      ) : (
        <Moon className="w-4 h-4" />
      )}
    </Button>
  )
}

export default SettingsPanel
