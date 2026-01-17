import { create } from 'zustand'

/**
 * [Dev Note: State Management]
 * ADC Builder 상태를 전역으로 관리하여 Step 간 이동 시 데이터 유지
 * Step 2에서 Step 1으로 뒤로 가기해도 입력값이 유지됨
 */

// 유효한 Amino Acid 문자 (FASTA 형식)
const VALID_AA_CHARS = /^[ACDEFGHIKLMNPQRSTVWY\s\n\r>]*$/i

export interface ADCBuilderFormData {
    // Step 1: Target & Antibody
    antibodyType: string
    customSequence: string
    targetName: string

    // Step 2: Payload & Linker
    payloadId: string
    linkerId: string
    dar: number

    // Step 3: Configuration
    mode: 'fast' | 'deep'
    jobName: string

    // Validation State
    sequenceError: string | null
}

interface ADCBuilderStore extends ADCBuilderFormData {
    // Actions
    setField: <K extends keyof ADCBuilderFormData>(key: K, value: ADCBuilderFormData[K]) => void
    validateSequence: (sequence: string) => boolean
    resetForm: () => void
}

const initialState: ADCBuilderFormData = {
    antibodyType: '',
    customSequence: '',
    targetName: '',
    payloadId: '',
    linkerId: '',
    dar: 4,
    mode: 'deep',
    jobName: '',
    sequenceError: null,
}

export const useADCBuilderStore = create<ADCBuilderStore>((set) => ({
    ...initialState,

    setField: (key, value) => set({ [key]: value } as Partial<ADCBuilderFormData>),

    /**
     * FASTA 서열 유효성 검사 (실시간)
     * - 유효하지 않은 문자가 있으면 즉시 경고
     * - 30분 돌리고 "서열 오류" 방지
     */
    validateSequence: (sequence: string) => {
        // 빈 문자열은 유효
        if (!sequence.trim()) {
            set({ sequenceError: null, customSequence: sequence })
            return true
        }

        // FASTA 헤더 제거 후 검사
        const lines = sequence.split('\n')
        const sequenceOnly = lines
            .filter(line => !line.startsWith('>'))
            .join('')
            .replace(/\s/g, '')

        if (!VALID_AA_CHARS.test(sequenceOnly)) {
            // 유효하지 않은 문자 찾기
            const invalidChars = sequenceOnly
                .split('')
                .filter(char => !/[ACDEFGHIKLMNPQRSTVWY]/i.test(char))
                .filter((v, i, a) => a.indexOf(v) === i) // unique
                .join(', ')

            set({
                sequenceError: `유효하지 않은 Amino Acid 코드: ${invalidChars}`,
                customSequence: sequence,
            })
            return false
        }

        // 최소 길이 검사
        if (sequenceOnly.length < 10) {
            set({
                sequenceError: '서열이 너무 짧습니다 (최소 10자)',
                customSequence: sequence,
            })
            return false
        }

        // 유효
        set({ sequenceError: null, customSequence: sequence })
        return true
    },

    resetForm: () => set(initialState),
}))
