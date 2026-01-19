export function TrustedBySection() {
    const logos = [
        "Pfizer", "Roche", "Novartis", "Merck", "Sanofi", "GSK", "AbbVie", "AstraZeneca"
    ]

    return (
        <section className="py-12 border-y border-white/5 bg-[#0F172A]/50 overflow-hidden">
            <div className="max-w-7xl mx-auto px-4 mb-8 text-center">
                <p className="text-sm font-semibold text-gray-500 uppercase tracking-widest">
                    Trusted Data Sources & Partners
                </p>
            </div>

            <div className="relative flex overflow-x-hidden group">
                <div className="animate-marquee whitespace-nowrap flex gap-16 items-center">
                    {logos.map((logo, index) => (
                        <span key={index} className="text-2xl font-bold text-gray-600 hover:text-gray-400 transition-colors cursor-default">
                            {logo}
                        </span>
                    ))}
                    {logos.map((logo, index) => (
                        <span key={`duplicate-${index}`} className="text-2xl font-bold text-gray-600 hover:text-gray-400 transition-colors cursor-default">
                            {logo}
                        </span>
                    ))}
                    {logos.map((logo, index) => (
                        <span key={`duplicate-2-${index}`} className="text-2xl font-bold text-gray-600 hover:text-gray-400 transition-colors cursor-default">
                            {logo}
                        </span>
                    ))}
                </div>
            </div>
        </section>
    )
}
