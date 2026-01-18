import React from 'react'

interface ContainerProps {
    children: React.ReactNode
    className?: string
    style?: React.CSSProperties
}

export function Container({ children, className = '', style = {} }: ContainerProps) {
    return (
        <div
            className={`px-4 sm:px-6 lg:px-8 ${className}`}
            style={{
                width: '100%',
                maxWidth: '80rem', // max-w-7xl equivalent
                marginLeft: 'auto',
                marginRight: 'auto',
                ...style
            }}
        >
            {children}
        </div>
    )
}
