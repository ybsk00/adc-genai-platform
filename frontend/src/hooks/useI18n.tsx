/**
 * useI18n Hook
 * React hook for translations (English-only)
 */
import { createContext, useContext } from 'react'
import type { ReactNode } from 'react'
import { getTranslations, t as translate } from '@/lib/i18n'
import type { Locale, Translations } from '@/lib/i18n'

interface I18nContextValue {
  locale: Locale
  t: Translations
  translate: (key: string) => string
}

const I18nContext = createContext<I18nContextValue | null>(null)

/**
 * I18n Provider
 * Provides translations context to the app
 */
export function I18nProvider({ children }: { children: ReactNode }) {
  const value: I18nContextValue = {
    locale: 'en',
    t: getTranslations(),
    translate: (key: string) => translate(key)
  }

  return (
    <I18nContext.Provider value={value}>
      {children}
    </I18nContext.Provider>
  )
}

/**
 * useI18n Hook
 * Get translations in components
 */
export function useI18n(): I18nContextValue {
  const context = useContext(I18nContext)

  if (!context) {
    // Provide default value if used outside provider
    return {
      locale: 'en',
      t: getTranslations(),
      translate: (key: string) => translate(key)
    }
  }

  return context
}

/**
 * useLocale Hook
 * Get current locale (always 'en')
 */
export function useLocale(): Locale {
  return 'en'
}

export default useI18n
