/**
 * Smart Disease Input Component
 * One-Click ADC Navigator - Disease → Target Selection
 *
 * 1단계: 질환명 검색/자동완성
 * 2단계: 타겟 선택 (HER2, TROP2, HER3 등 golden_set 기반)
 */

import { useState, useEffect, useCallback, useRef } from 'react';
import { Search, Target, FlaskConical, Loader2, Rocket, X, ChevronRight, ArrowLeft } from 'lucide-react';
import { Input } from '@/components/ui/input';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { cn } from '@/lib/utils';
import { motion, AnimatePresence } from 'framer-motion';
import { API_BASE_URL } from '@/lib/api';

interface DiseaseSuggestion {
  disease: string;
  primary_targets: string[];
  approved_adc_count: number;
  avg_orr: number | null;
}

interface SmartDiseaseInputProps {
  onSubmit: (diseaseName: string, targetProtein: string) => void;
  isLoading?: boolean;
  disabled?: boolean;
}

// Debounce hook
function useDebounce<T>(value: T, delay: number): T {
  const [debouncedValue, setDebouncedValue] = useState<T>(value);

  useEffect(() => {
    const handler = setTimeout(() => {
      setDebouncedValue(value);
    }, delay);

    return () => clearTimeout(handler);
  }, [value, delay]);

  return debouncedValue;
}

// 타겟별 아이콘 색상 매핑
const TARGET_COLORS: Record<string, string> = {
  HER2: 'from-pink-500 to-rose-500',
  TROP2: 'from-emerald-500 to-teal-500',
  HER3: 'from-purple-500 to-violet-500',
  'Nectin-4': 'from-blue-500 to-cyan-500',
  CD30: 'from-orange-500 to-amber-500',
  CD79b: 'from-indigo-500 to-blue-500',
  'c-Met': 'from-lime-500 to-green-500',
  'Tissue Factor': 'from-red-500 to-pink-500',
  CD33: 'from-yellow-500 to-orange-500',
  FRa: 'from-teal-500 to-cyan-500',
  gpNMB: 'from-gray-500 to-slate-500',
};

function getTargetColor(target: string): string {
  return TARGET_COLORS[target] || 'from-violet-500 to-indigo-500';
}

export function SmartDiseaseInput({
  onSubmit,
  isLoading = false,
  disabled = false,
}: SmartDiseaseInputProps) {
  const [query, setQuery] = useState('');
  const [suggestions, setSuggestions] = useState<DiseaseSuggestion[]>([]);
  const [isFetching, setIsFetching] = useState(false);
  const [showDropdown, setShowDropdown] = useState(false);
  const [selectedIndex, setSelectedIndex] = useState(-1);

  // 2단계: 타겟 선택
  const [selectedDisease, setSelectedDisease] = useState<DiseaseSuggestion | null>(null);

  const inputRef = useRef<HTMLInputElement>(null);
  const dropdownRef = useRef<HTMLDivElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);

  const debouncedQuery = useDebounce(query, 300);

  // Fetch suggestions
  const fetchSuggestions = useCallback(async (searchQuery: string) => {
    if (searchQuery.length < 2) {
      setSuggestions([]);
      return;
    }

    setIsFetching(true);
    try {
      const response = await fetch(
        `${API_BASE_URL}/api/design/navigator/suggest-disease?query=${encodeURIComponent(searchQuery)}&limit=8`
      );

      if (response.ok) {
        const data = await response.json();
        setSuggestions(data.suggestions || []);
        setShowDropdown(true);
      }
    } catch (error) {
      console.error('Failed to fetch disease suggestions:', error);
      setSuggestions([]);
    } finally {
      setIsFetching(false);
    }
  }, []);

  // Fetch when debounced query changes
  useEffect(() => {
    if (!selectedDisease) {
      fetchSuggestions(debouncedQuery);
    }
  }, [debouncedQuery, fetchSuggestions, selectedDisease]);

  // Handle outside click
  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (
        containerRef.current &&
        !containerRef.current.contains(event.target as Node)
      ) {
        setShowDropdown(false);
      }
    };

    document.addEventListener('mousedown', handleClickOutside);
    return () => document.removeEventListener('mousedown', handleClickOutside);
  }, []);

  // Handle keyboard navigation
  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (selectedDisease) return; // 타겟 선택 단계에서는 키보드 무시

    if (!showDropdown || suggestions.length === 0) {
      if (e.key === 'Enter' && query.trim()) {
        e.preventDefault();
      }
      return;
    }

    switch (e.key) {
      case 'ArrowDown':
        e.preventDefault();
        setSelectedIndex((prev) =>
          prev < suggestions.length - 1 ? prev + 1 : 0
        );
        break;
      case 'ArrowUp':
        e.preventDefault();
        setSelectedIndex((prev) =>
          prev > 0 ? prev - 1 : suggestions.length - 1
        );
        break;
      case 'Enter':
        e.preventDefault();
        if (selectedIndex >= 0) {
          handleSelectDisease(suggestions[selectedIndex]);
        }
        break;
      case 'Escape':
        setShowDropdown(false);
        setSelectedIndex(-1);
        break;
    }
  };

  // 1단계: 질환 선택 → 타겟 선택 UI 표시
  const handleSelectDisease = (suggestion: DiseaseSuggestion) => {
    setQuery(suggestion.disease);
    setShowDropdown(false);
    setSelectedIndex(-1);
    setSelectedDisease(suggestion);
  };

  // 2단계: 타겟 선택 → Navigator 실행
  const handleSelectTarget = (target: string) => {
    if (selectedDisease && !disabled && !isLoading) {
      onSubmit(selectedDisease.disease, target);
    }
  };

  // 뒤로가기
  const handleBackToSearch = () => {
    setSelectedDisease(null);
    setQuery('');
    setSuggestions([]);
    inputRef.current?.focus();
  };

  const handleClear = () => {
    setQuery('');
    setSuggestions([]);
    setShowDropdown(false);
    setSelectedDisease(null);
    inputRef.current?.focus();
  };

  return (
    <div ref={containerRef} className="relative w-full max-w-2xl mx-auto">
      <AnimatePresence mode="wait">
        {!selectedDisease ? (
          /* ============ 1단계: 질환 검색 ============ */
          <motion.div
            key="search"
            initial={{ opacity: 0, y: -10 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -10 }}
          >
            {/* Main Input Container */}
            <div className="relative">
              <div className="flex items-center gap-2 p-2 bg-gradient-to-r from-slate-900 to-slate-800 rounded-xl border border-slate-700 shadow-lg">
                <div className="pl-3 flex items-center">
                  {isFetching ? (
                    <Loader2 className="h-5 w-5 text-slate-400 animate-spin" />
                  ) : (
                    <Search className="h-5 w-5 text-slate-400" />
                  )}
                </div>

                <Input
                  ref={inputRef}
                  type="text"
                  value={query}
                  onChange={(e) => setQuery(e.target.value)}
                  onKeyDown={handleKeyDown}
                  onFocus={() => suggestions.length > 0 && setShowDropdown(true)}
                  placeholder="Enter disease name (e.g., Breast Cancer, Lymphoma)"
                  className="flex-1 bg-transparent border-none text-white placeholder:text-slate-500 focus-visible:ring-0 focus-visible:ring-offset-0 text-lg"
                  disabled={disabled || isLoading}
                />

                {query && (
                  <Button
                    variant="ghost"
                    size="icon"
                    onClick={handleClear}
                    className="h-8 w-8 text-slate-400 hover:text-white"
                  >
                    <X className="h-4 w-4" />
                  </Button>
                )}
              </div>
            </div>

            {/* Suggestions Dropdown */}
            {showDropdown && suggestions.length > 0 && (
              <div
                ref={dropdownRef}
                className="absolute z-50 w-full mt-2 bg-slate-900 border border-slate-700 rounded-xl shadow-2xl overflow-hidden"
              >
                <div className="p-2 border-b border-slate-700">
                  <span className="text-xs text-slate-500 font-medium">
                    {suggestions.length} disease{suggestions.length > 1 ? 's' : ''} found — click to select target
                  </span>
                </div>
                <ul className="max-h-80 overflow-y-auto">
                  {suggestions.map((suggestion, index) => (
                    <li
                      key={suggestion.disease}
                      onClick={() => handleSelectDisease(suggestion)}
                      className={cn(
                        'px-4 py-3 cursor-pointer transition-colors border-b border-slate-800 last:border-0',
                        selectedIndex === index
                          ? 'bg-slate-800'
                          : 'hover:bg-slate-800/50'
                      )}
                    >
                      <div className="flex items-center justify-between gap-4">
                        <div className="flex-1">
                          <div className="flex items-center gap-2">
                            <span className="text-white font-medium">
                              {suggestion.disease}
                            </span>
                            {suggestion.approved_adc_count > 0 && (
                              <Badge
                                variant="secondary"
                                className="bg-emerald-500/20 text-emerald-400 text-xs"
                              >
                                {suggestion.approved_adc_count} Approved ADC{suggestion.approved_adc_count > 1 ? 's' : ''}
                              </Badge>
                            )}
                          </div>

                          {suggestion.primary_targets.length > 0 && (
                            <div className="flex items-center gap-2 mt-1.5">
                              <Target className="h-3.5 w-3.5 text-slate-500" />
                              <span className="text-sm text-slate-400">
                                Targets:{' '}
                                <span className="text-violet-400">
                                  {suggestion.primary_targets.join(', ')}
                                </span>
                              </span>
                            </div>
                          )}
                        </div>

                        <ChevronRight className="h-5 w-5 text-slate-500" />
                      </div>
                    </li>
                  ))}
                </ul>
              </div>
            )}

            {/* Helper Text */}
            <p className="text-center text-sm text-slate-500 mt-3">
              <FlaskConical className="inline h-4 w-4 mr-1" />
              Search a disease, then select a molecular target from approved Golden Set ADCs
            </p>
          </motion.div>
        ) : (
          /* ============ 2단계: 타겟 선택 ============ */
          <motion.div
            key="target-select"
            initial={{ opacity: 0, y: 10 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: 10 }}
          >
            {/* Header: 선택된 질환 */}
            <div className="flex items-center gap-3 mb-6">
              <Button
                variant="ghost"
                size="icon"
                onClick={handleBackToSearch}
                className="h-8 w-8 text-slate-400 hover:text-white"
                disabled={isLoading}
              >
                <ArrowLeft className="h-4 w-4" />
              </Button>
              <div className="flex-1">
                <p className="text-sm text-slate-500">Selected Disease</p>
                <div className="flex items-center gap-2">
                  <h3 className="text-xl font-bold text-white">{selectedDisease.disease}</h3>
                  {selectedDisease.approved_adc_count > 0 && (
                    <Badge className="bg-emerald-500/20 text-emerald-400 text-xs">
                      {selectedDisease.approved_adc_count} Approved ADCs
                    </Badge>
                  )}
                </div>
              </div>
              {selectedDisease.avg_orr !== null && (
                <div className="text-right">
                  <p className="text-xs text-slate-500">Avg ORR</p>
                  <p className="text-lg font-bold text-amber-400">
                    {Number(selectedDisease.avg_orr ?? 0).toFixed(1)}%
                  </p>
                </div>
              )}
            </div>

            {/* 타겟 카드 그리드 */}
            <div className="mb-4">
              <p className="text-sm text-slate-400 mb-3 flex items-center gap-2">
                <Target className="h-4 w-4 text-violet-400" />
                Select a molecular target to design your ADC
              </p>

              <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-3">
                {selectedDisease.primary_targets.map((target, index) => (
                  <motion.button
                    key={target}
                    initial={{ opacity: 0, scale: 0.9 }}
                    animate={{ opacity: 1, scale: 1 }}
                    transition={{ delay: index * 0.08 }}
                    onClick={() => handleSelectTarget(target)}
                    disabled={disabled || isLoading}
                    className={cn(
                      'group relative p-5 rounded-xl border-2 text-left transition-all',
                      'bg-slate-900/80 border-slate-700',
                      'hover:border-violet-500/70 hover:bg-slate-800/80 hover:shadow-lg hover:shadow-violet-500/10',
                      'focus:outline-none focus:ring-2 focus:ring-violet-500/50',
                      'disabled:opacity-50 disabled:cursor-not-allowed'
                    )}
                  >
                    {/* Target Icon */}
                    <div className={cn(
                      'w-10 h-10 rounded-lg bg-gradient-to-r flex items-center justify-center mb-3',
                      getTargetColor(target)
                    )}>
                      <Target className="h-5 w-5 text-white" />
                    </div>

                    {/* Target Name */}
                    <h4 className="text-lg font-bold text-white group-hover:text-violet-300 transition-colors">
                      {target}
                    </h4>

                    {/* Arrow */}
                    <div className="absolute top-5 right-4 opacity-0 group-hover:opacity-100 transition-opacity">
                      <Rocket className="h-5 w-5 text-violet-400" />
                    </div>
                  </motion.button>
                ))}
              </div>
            </div>

            {/* Loading State */}
            {isLoading && (
              <div className="flex items-center justify-center gap-2 text-violet-400 mt-4">
                <Loader2 className="h-5 w-5 animate-spin" />
                <span>Starting Navigator...</span>
              </div>
            )}
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  );
}

export default SmartDiseaseInput;
