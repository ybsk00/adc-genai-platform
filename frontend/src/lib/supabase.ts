import { createClient, SupabaseClient } from '@supabase/supabase-js'

const supabaseUrl = import.meta.env.VITE_SUPABASE_URL
const supabaseAnonKey = import.meta.env.VITE_SUPABASE_ANON_KEY

// Supabase 클라이언트를 조건부로 생성 (환경변수가 없으면 null)
let supabaseClient: SupabaseClient | null = null

if (supabaseUrl && supabaseAnonKey) {
    supabaseClient = createClient(supabaseUrl, supabaseAnonKey)
} else {
    console.warn('Supabase environment variables are not set. Auth features will not work.')
}

export const supabase = supabaseClient

// Auth helper functions
export const signInWithEmail = async (email: string, password: string) => {
    if (!supabase) {
        return { data: null, error: { message: 'Supabase is not configured' } }
    }
    const { data, error } = await supabase.auth.signInWithPassword({
        email,
        password,
    })
    return { data, error }
}

export const signUpWithEmail = async (email: string, password: string, fullName?: string) => {
    if (!supabase) {
        return { data: null, error: { message: 'Supabase is not configured' } }
    }
    const { data, error } = await supabase.auth.signUp({
        email,
        password,
        options: {
            data: {
                full_name: fullName,
            },
        },
    })
    return { data, error }
}

export const signInWithGoogle = async () => {
    if (!supabase) {
        return { data: null, error: { message: 'Supabase is not configured' } }
    }
    const { data, error } = await supabase.auth.signInWithOAuth({
        provider: 'google',
        options: {
            redirectTo: `${window.location.origin}/dashboard`,
        },
    })
    return { data, error }
}

export const signOut = async () => {
    if (!supabase) {
        return { error: { message: 'Supabase is not configured' } }
    }
    const { error } = await supabase.auth.signOut()
    return { error }
}

export const getSession = async () => {
    if (!supabase) {
        return { session: null, error: { message: 'Supabase is not configured' } }
    }
    const { data: { session }, error } = await supabase.auth.getSession()
    return { session, error }
}

export const getUser = async () => {
    if (!supabase) {
        return { user: null, error: { message: 'Supabase is not configured' } }
    }
    const { data: { user }, error } = await supabase.auth.getUser()
    return { user, error }
}

