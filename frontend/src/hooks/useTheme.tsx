/**
 * useTheme Hook
 * 다크 모드 / 라이트 모드 테마 관리
 */
import { useState, useEffect, useCallback, createContext, useContext } from 'react'
import type { ReactNode } from 'react'

export type Theme = 'light' | 'dark' | 'system'

interface ThemeContextValue {
  theme: Theme
  resolvedTheme: 'light' | 'dark'
  setTheme: (theme: Theme) => void
  toggleTheme: () => void
}

const THEME_STORAGE_KEY = 'adc-platform-theme'

const ThemeContext = createContext<ThemeContextValue | null>(null)

/**
 * 시스템 테마 감지
 */
function getSystemTheme(): 'light' | 'dark' {
  if (typeof window === 'undefined') return 'light'
  return window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light'
}

/**
 * 저장된 테마 가져오기
 */
function getStoredTheme(): Theme {
  if (typeof window === 'undefined') return 'system'
  const stored = localStorage.getItem(THEME_STORAGE_KEY)
  if (stored === 'light' || stored === 'dark' || stored === 'system') {
    return stored
  }
  return 'system'
}

/**
 * Theme Provider
 */
export function ThemeProvider({ children }: { children: ReactNode }) {
  const [theme, setThemeState] = useState<Theme>(getStoredTheme)
  const [resolvedTheme, setResolvedTheme] = useState<'light' | 'dark'>(
    theme === 'system' ? getSystemTheme() : theme
  )

  // 테마 변경
  const setTheme = useCallback((newTheme: Theme) => {
    localStorage.setItem(THEME_STORAGE_KEY, newTheme)
    setThemeState(newTheme)

    const resolved = newTheme === 'system' ? getSystemTheme() : newTheme
    setResolvedTheme(resolved)

    // DOM 업데이트
    updateDOM(resolved)
  }, [])

  // 테마 토글
  const toggleTheme = useCallback(() => {
    const newTheme = resolvedTheme === 'dark' ? 'light' : 'dark'
    setTheme(newTheme)
  }, [resolvedTheme, setTheme])

  // 시스템 테마 변경 감지
  useEffect(() => {
    const mediaQuery = window.matchMedia('(prefers-color-scheme: dark)')

    const handleChange = (e: MediaQueryListEvent) => {
      if (theme === 'system') {
        const newResolved = e.matches ? 'dark' : 'light'
        setResolvedTheme(newResolved)
        updateDOM(newResolved)
      }
    }

    mediaQuery.addEventListener('change', handleChange)
    return () => mediaQuery.removeEventListener('change', handleChange)
  }, [theme])

  // 초기 DOM 업데이트
  useEffect(() => {
    updateDOM(resolvedTheme)
  }, [])

  const value: ThemeContextValue = {
    theme,
    resolvedTheme,
    setTheme,
    toggleTheme
  }

  return (
    <ThemeContext.Provider value={value}>
      {children}
    </ThemeContext.Provider>
  )
}

/**
 * DOM 클래스 업데이트
 */
function updateDOM(theme: 'light' | 'dark') {
  const root = document.documentElement

  if (theme === 'dark') {
    root.classList.add('dark')
    root.style.colorScheme = 'dark'
  } else {
    root.classList.remove('dark')
    root.style.colorScheme = 'light'
  }
}

/**
 * useTheme Hook
 */
export function useTheme(): ThemeContextValue {
  const context = useContext(ThemeContext)

  if (!context) {
    // Context 없이 사용 시 기본값 제공
    const theme = getStoredTheme()
    const resolved = theme === 'system' ? getSystemTheme() : theme

    return {
      theme,
      resolvedTheme: resolved,
      setTheme: (newTheme: Theme) => {
        localStorage.setItem(THEME_STORAGE_KEY, newTheme)
        updateDOM(newTheme === 'system' ? getSystemTheme() : newTheme)
      },
      toggleTheme: () => {
        const newTheme = resolved === 'dark' ? 'light' : 'dark'
        localStorage.setItem(THEME_STORAGE_KEY, newTheme)
        updateDOM(newTheme)
      }
    }
  }

  return context
}

export default useTheme
