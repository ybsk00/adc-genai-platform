/**
 * Smart Disease Input Component
 * One-Click ADC Navigator - Disease Name Autocomplete
 *
 * 질환명 자동완성 및 타겟 제안 컴포넌트
 */

import { useState, useEffect, useCallback, useRef } from 'react';
import { Search, Target, FlaskConical, Loader2, Rocket, X } from 'lucide-react';
import { Input } from '@/components/ui/input';
import { Button } from '@/components/ui/button';
import { Badge } from '@/components/ui/badge';
import { cn } from '@/lib/utils';
import { API_BASE_URL } from '@/lib/api';

interface DiseaseSuggestion {
  disease: string;
  primary_targets: string[];
  approved_adc_count: number;
  avg_orr: number | null;
}

interface SmartDiseaseInputProps {
  onSelect?: (suggestion: DiseaseSuggestion) => void;
  onSubmit: (diseaseName: string) => void;
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

export function SmartDiseaseInput({
  onSelect,
  onSubmit,
  isLoading = false,
  disabled = false,
}: SmartDiseaseInputProps) {
  const [query, setQuery] = useState('');
  const [suggestions, setSuggestions] = useState<DiseaseSuggestion[]>([]);
  const [isFetching, setIsFetching] = useState(false);
  const [showDropdown, setShowDropdown] = useState(false);
  const [selectedIndex, setSelectedIndex] = useState(-1);

  const inputRef = useRef<HTMLInputElement>(null);
  const dropdownRef = useRef<HTMLDivElement>(null);

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
    fetchSuggestions(debouncedQuery);
  }, [debouncedQuery, fetchSuggestions]);

  // Handle outside click
  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (
        dropdownRef.current &&
        !dropdownRef.current.contains(event.target as Node) &&
        inputRef.current &&
        !inputRef.current.contains(event.target as Node)
      ) {
        setShowDropdown(false);
      }
    };

    document.addEventListener('mousedown', handleClickOutside);
    return () => document.removeEventListener('mousedown', handleClickOutside);
  }, []);

  // Handle keyboard navigation
  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (!showDropdown || suggestions.length === 0) {
      if (e.key === 'Enter' && query.trim()) {
        e.preventDefault();
        handleSubmit();
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
          handleSelectSuggestion(suggestions[selectedIndex]);
        } else if (query.trim()) {
          handleSubmit();
        }
        break;
      case 'Escape':
        setShowDropdown(false);
        setSelectedIndex(-1);
        break;
    }
  };

  const handleSelectSuggestion = (suggestion: DiseaseSuggestion) => {
    setQuery(suggestion.disease);
    setShowDropdown(false);
    setSelectedIndex(-1);
    onSelect?.(suggestion);
  };

  const handleSubmit = () => {
    if (query.trim() && !disabled && !isLoading) {
      onSubmit(query.trim());
    }
  };

  const handleClear = () => {
    setQuery('');
    setSuggestions([]);
    setShowDropdown(false);
    inputRef.current?.focus();
  };

  return (
    <div className="relative w-full max-w-2xl mx-auto">
      {/* Main Input Container */}
      <div className="relative">
        <div className="flex items-center gap-2 p-2 bg-gradient-to-r from-slate-900 to-slate-800 rounded-xl border border-slate-700 shadow-lg">
          {/* Search Icon */}
          <div className="pl-3 flex items-center">
            {isFetching ? (
              <Loader2 className="h-5 w-5 text-slate-400 animate-spin" />
            ) : (
              <Search className="h-5 w-5 text-slate-400" />
            )}
          </div>

          {/* Input */}
          <Input
            ref={inputRef}
            type="text"
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            onKeyDown={handleKeyDown}
            onFocus={() => suggestions.length > 0 && setShowDropdown(true)}
            placeholder="Enter disease name (e.g., Pancreatic Cancer, HER2+ Breast Cancer)"
            className="flex-1 bg-transparent border-none text-white placeholder:text-slate-500 focus-visible:ring-0 focus-visible:ring-offset-0 text-lg"
            disabled={disabled || isLoading}
          />

          {/* Clear Button */}
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

          {/* Submit Button */}
          <Button
            onClick={handleSubmit}
            disabled={!query.trim() || disabled || isLoading}
            className="bg-gradient-to-r from-violet-600 to-indigo-600 hover:from-violet-700 hover:to-indigo-700 text-white px-6"
          >
            {isLoading ? (
              <>
                <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                Processing...
              </>
            ) : (
              <>
                <Rocket className="mr-2 h-4 w-4" />
                One-Click Navigate
              </>
            )}
          </Button>
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
              {suggestions.length} disease{suggestions.length > 1 ? 's' : ''} found
            </span>
          </div>
          <ul className="max-h-80 overflow-y-auto">
            {suggestions.map((suggestion, index) => (
              <li
                key={suggestion.disease}
                onClick={() => handleSelectSuggestion(suggestion)}
                className={cn(
                  'px-4 py-3 cursor-pointer transition-colors border-b border-slate-800 last:border-0',
                  selectedIndex === index
                    ? 'bg-slate-800'
                    : 'hover:bg-slate-800/50'
                )}
              >
                <div className="flex items-start justify-between gap-4">
                  {/* Disease Name */}
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

                    {/* Primary Targets */}
                    {suggestion.primary_targets.length > 0 && (
                      <div className="flex items-center gap-2 mt-1.5">
                        <Target className="h-3.5 w-3.5 text-slate-500" />
                        <span className="text-sm text-slate-400">
                          Targets:{' '}
                          <span className="text-violet-400">
                            {suggestion.primary_targets.slice(0, 3).join(', ')}
                            {suggestion.primary_targets.length > 3 && ' ...'}
                          </span>
                        </span>
                      </div>
                    )}
                  </div>

                  {/* ORR Badge */}
                  {suggestion.avg_orr !== null && (
                    <div className="flex flex-col items-end">
                      <span className="text-xs text-slate-500">Avg ORR</span>
                      <span className="text-sm font-semibold text-amber-400">
                        {suggestion.avg_orr.toFixed(1)}%
                      </span>
                    </div>
                  )}
                </div>
              </li>
            ))}
          </ul>
        </div>
      )}

      {/* Helper Text */}
      <p className="text-center text-sm text-slate-500 mt-3">
        <FlaskConical className="inline h-4 w-4 mr-1" />
        Our AI will analyze 7,600+ compounds to design the optimal ADC for your target disease
      </p>
    </div>
  );
}

export default SmartDiseaseInput;
