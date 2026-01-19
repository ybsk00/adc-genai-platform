import { useState } from 'react'
import { motion } from 'framer-motion'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Badge } from '@/components/ui/badge'
import {
    Select,
    SelectContent,
    SelectItem,
    SelectTrigger,
    SelectValue,
} from '@/components/ui/select'
import {
    Dialog,
    DialogContent,
    DialogDescription,
    DialogHeader,
    DialogTitle,
} from '@/components/ui/dialog'
import {
    Search,
    Filter,
    FlaskConical,
    Building2,
    Calendar,
    Target,
    Pill,
    Link2,
    ExternalLink,
    ChevronRight
} from 'lucide-react'

// FDA 승인 ADC 데이터 (Golden Set)
const goldenSetData = [
    {
        id: 'gs_1',
        drugName: 'Kadcyla',
        brandName: 'Kadcyla (ado-trastuzumab emtansine)',
        target: 'HER2',
        antibody: 'Trastuzumab',
        payload: 'DM1',
        linker: 'MCC (non-cleavable)',
        dar: 3.5,
        company: 'Roche/Genentech',
        approvalDate: '2013-02-22',
        indications: ['HER2+ metastatic breast cancer'],
        status: 'Approved'
    },
    {
        id: 'gs_2',
        drugName: 'Enhertu',
        brandName: 'Enhertu (trastuzumab deruxtecan)',
        target: 'HER2',
        antibody: 'Trastuzumab',
        payload: 'DXd',
        linker: 'GGFG (cleavable)',
        dar: 8.0,
        company: 'Daiichi Sankyo/AstraZeneca',
        approvalDate: '2019-12-20',
        indications: ['HER2+ breast cancer', 'HER2-low breast cancer', 'HER2+ gastric cancer'],
        status: 'Approved'
    },
    {
        id: 'gs_3',
        drugName: 'Trodelvy',
        brandName: 'Trodelvy (sacituzumab govitecan)',
        target: 'TROP-2',
        antibody: 'Sacituzumab',
        payload: 'SN-38',
        linker: 'CL2A (cleavable)',
        dar: 7.6,
        company: 'Gilead Sciences',
        approvalDate: '2020-04-22',
        indications: ['Triple-negative breast cancer', 'HR+/HER2- breast cancer'],
        status: 'Approved'
    },
    {
        id: 'gs_4',
        drugName: 'Adcetris',
        brandName: 'Adcetris (brentuximab vedotin)',
        target: 'CD30',
        antibody: 'Brentuximab',
        payload: 'MMAE',
        linker: 'Val-Cit (cleavable)',
        dar: 4.0,
        company: 'Seagen',
        approvalDate: '2011-08-19',
        indications: ['Hodgkin lymphoma', 'ALCL'],
        status: 'Approved'
    },
    {
        id: 'gs_5',
        drugName: 'Padcev',
        brandName: 'Padcev (enfortumab vedotin)',
        target: 'Nectin-4',
        antibody: 'Enfortumab',
        payload: 'MMAE',
        linker: 'Val-Cit (cleavable)',
        dar: 3.8,
        company: 'Seagen/Astellas',
        approvalDate: '2019-12-18',
        indications: ['Urothelial cancer'],
        status: 'Approved'
    },
    {
        id: 'gs_6',
        drugName: 'Polivy',
        brandName: 'Polivy (polatuzumab vedotin)',
        target: 'CD79b',
        antibody: 'Polatuzumab',
        payload: 'MMAE',
        linker: 'Val-Cit (cleavable)',
        dar: 3.5,
        company: 'Roche',
        approvalDate: '2019-06-10',
        indications: ['DLBCL'],
        status: 'Approved'
    },
]

// 고유 타겟 목록
const targets = [...new Set(goldenSetData.map(d => d.target))]
const payloads = [...new Set(goldenSetData.map(d => d.payload))]

interface GoldenSetItem {
    id: string
    drugName: string
    brandName: string
    target: string
    antibody: string
    payload: string
    linker: string
    dar: number
    company: string
    approvalDate: string
    indications: string[]
    status: string
}

export function GoldenSetLibrary() {
    const [searchQuery, setSearchQuery] = useState('')
    const [targetFilter, setTargetFilter] = useState<string>('all')
    const [payloadFilter, setPayloadFilter] = useState<string>('all')
    const [selectedItem, setSelectedItem] = useState<GoldenSetItem | null>(null)
    const [detailOpen, setDetailOpen] = useState(false)

    // 필터링
    const filteredData = goldenSetData.filter(item => {
        const matchesSearch =
            item.drugName.toLowerCase().includes(searchQuery.toLowerCase()) ||
            item.brandName.toLowerCase().includes(searchQuery.toLowerCase()) ||
            item.company.toLowerCase().includes(searchQuery.toLowerCase())

        const matchesTarget = targetFilter === 'all' || item.target === targetFilter
        const matchesPayload = payloadFilter === 'all' || item.payload === payloadFilter

        return matchesSearch && matchesTarget && matchesPayload
    })

    const handleViewDetail = (item: GoldenSetItem) => {
        setSelectedItem(item)
        setDetailOpen(true)
    }

    const handleUseAsTemplate = (item: GoldenSetItem) => {
        // TODO: ADC Builder로 이동하면서 데이터 전달
        window.location.href = `/dashboard/builder?template=${item.id}`
    }

    return (
        <div className="space-y-6">
            {/* Page Header */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
            >
                <h1 className="text-2xl font-bold text-white">Golden Set Library</h1>
                <p className="text-slate-400 mt-1">FDA Approved ADC Database - 15 Verified ADC Entries</p>
            </motion.div>

            {/* Search & Filters */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.1 }}
            >
                <Card className="bg-slate-900 border-slate-800">
                    <CardContent className="p-4">
                        <div className="flex flex-col md:flex-row gap-4">
                            {/* Search */}
                            <div className="relative flex-1">
                                <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-slate-500" />
                                <Input
                                    placeholder="Search by drug name, company..."
                                    value={searchQuery}
                                    onChange={(e) => setSearchQuery(e.target.value)}
                                    className="pl-10 bg-slate-950 border-slate-800 text-white placeholder:text-slate-500"
                                />
                            </div>

                            {/* Target Filter */}
                            <Select value={targetFilter} onValueChange={setTargetFilter}>
                                <SelectTrigger className="w-full md:w-[180px] bg-slate-950 border-slate-800 text-white">
                                    <Target className="w-4 h-4 mr-2 text-slate-400" />
                                    <SelectValue placeholder="Target" />
                                </SelectTrigger>
                                <SelectContent className="bg-slate-900 border-slate-800 text-white">
                                    <SelectItem value="all" className="focus:bg-slate-800 focus:text-white">All Targets</SelectItem>
                                    {targets.map(t => (
                                        <SelectItem key={t} value={t} className="focus:bg-slate-800 focus:text-white">{t}</SelectItem>
                                    ))}
                                </SelectContent>
                            </Select>

                            {/* Payload Filter */}
                            <Select value={payloadFilter} onValueChange={setPayloadFilter}>
                                <SelectTrigger className="w-full md:w-[180px] bg-slate-950 border-slate-800 text-white">
                                    <Pill className="w-4 h-4 mr-2 text-slate-400" />
                                    <SelectValue placeholder="Payload" />
                                </SelectTrigger>
                                <SelectContent className="bg-slate-900 border-slate-800 text-white">
                                    <SelectItem value="all" className="focus:bg-slate-800 focus:text-white">All Payloads</SelectItem>
                                    {payloads.map(p => (
                                        <SelectItem key={p} value={p} className="focus:bg-slate-800 focus:text-white">{p}</SelectItem>
                                    ))}
                                </SelectContent>
                            </Select>
                        </div>
                    </CardContent>
                </Card>
            </motion.div>

            {/* Results Count */}
            <div className="flex items-center justify-between">
                <p className="text-sm text-slate-500">
                    {filteredData.length} results
                </p>
                <Button variant="outline" size="sm" className="border-slate-700 text-slate-300 hover:bg-slate-800 hover:text-white">
                    <Filter className="w-4 h-4 mr-2" />
                    Advanced Filter
                </Button>
            </div>

            {/* ADC Cards Grid */}
            <motion.div
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                transition={{ delay: 0.2 }}
                className="grid md:grid-cols-2 lg:grid-cols-3 gap-4"
            >
                {filteredData.map((item, index) => (
                    <motion.div
                        key={item.id}
                        initial={{ opacity: 0, y: 20 }}
                        animate={{ opacity: 1, y: 0 }}
                        transition={{ delay: 0.1 + index * 0.05 }}
                    >
                        <Card className="bg-slate-900 border-slate-800 hover:border-slate-700 hover:shadow-lg hover:shadow-blue-500/5 transition-all cursor-pointer group">
                            <CardHeader className="pb-3">
                                <div className="flex items-start justify-between">
                                    <div>
                                        <CardTitle className="text-lg text-white group-hover:text-blue-400 transition-colors">
                                            {item.drugName}
                                        </CardTitle>
                                        <CardDescription className="text-xs mt-1 text-slate-400">
                                            {item.brandName}
                                        </CardDescription>
                                    </div>
                                    <Badge className="bg-green-500/10 text-green-400 border-green-500/20 hover:bg-green-500/20">
                                        {item.status}
                                    </Badge>
                                </div>
                            </CardHeader>
                            <CardContent className="space-y-3">
                                <div className="grid grid-cols-2 gap-2 text-sm">
                                    <div className="flex items-center gap-2 text-slate-300">
                                        <Target className="w-4 h-4 text-blue-400" />
                                        <span>{item.target}</span>
                                    </div>
                                    <div className="flex items-center gap-2 text-slate-300">
                                        <Pill className="w-4 h-4 text-purple-400" />
                                        <span>{item.payload}</span>
                                    </div>
                                    <div className="flex items-center gap-2 text-slate-300">
                                        <Link2 className="w-4 h-4 text-orange-400" />
                                        <span className="truncate">{item.linker.split(' ')[0]}</span>
                                    </div>
                                    <div className="flex items-center gap-2 text-slate-300">
                                        <FlaskConical className="w-4 h-4 text-green-400" />
                                        <span>DAR {item.dar}</span>
                                    </div>
                                </div>

                                <div className="flex items-center gap-2 text-sm text-slate-500">
                                    <Building2 className="w-4 h-4" />
                                    <span className="truncate">{item.company}</span>
                                </div>

                                <div className="flex items-center gap-2 text-sm text-slate-500">
                                    <Calendar className="w-4 h-4" />
                                    <span>Approval Date: {item.approvalDate}</span>
                                </div>

                                <div className="flex gap-2 pt-2">
                                    <Button
                                        variant="outline"
                                        size="sm"
                                        className="flex-1 border-slate-700 text-slate-300 hover:bg-slate-800 hover:text-white"
                                        onClick={() => handleViewDetail(item)}
                                    >
                                        View Details
                                    </Button>
                                    <Button
                                        size="sm"
                                        className="flex-1 bg-blue-600 hover:bg-blue-700 text-white"
                                        onClick={() => handleUseAsTemplate(item)}
                                    >
                                        Use Template
                                        <ChevronRight className="w-4 h-4 ml-1" />
                                    </Button>
                                </div>
                            </CardContent>
                        </Card>
                    </motion.div>
                ))}
            </motion.div>

            {/* Detail Modal */}
            <Dialog open={detailOpen} onOpenChange={setDetailOpen}>
                <DialogContent className="max-w-2xl bg-slate-900 border-slate-800 text-white">
                    {selectedItem && (
                        <>
                            <DialogHeader>
                                <DialogTitle className="text-xl text-white">{selectedItem.drugName}</DialogTitle>
                                <DialogDescription className="text-slate-400">{selectedItem.brandName}</DialogDescription>
                            </DialogHeader>

                            <div className="space-y-6 py-4">
                                {/* Basic Info */}
                                <div className="grid grid-cols-2 gap-4">
                                    <div>
                                        <p className="text-sm text-slate-500">Target</p>
                                        <p className="font-medium text-slate-200">{selectedItem.target}</p>
                                    </div>
                                    <div>
                                        <p className="text-sm text-slate-500">Payload</p>
                                        <p className="font-medium text-slate-200">{selectedItem.payload}</p>
                                    </div>
                                    <div>
                                        <p className="text-sm text-slate-500">Antibody</p>
                                        <p className="font-medium text-slate-200">{selectedItem.antibody}</p>
                                    </div>
                                    <div>
                                        <p className="text-sm text-slate-500">Linker</p>
                                        <p className="font-medium text-slate-200">{selectedItem.linker}</p>
                                    </div>
                                    <div>
                                        <p className="text-sm text-slate-500">DAR</p>
                                        <p className="font-medium text-slate-200">{selectedItem.dar}</p>
                                    </div>
                                    <div>
                                        <p className="text-sm text-slate-500">Developer</p>
                                        <p className="font-medium text-slate-200">{selectedItem.company}</p>
                                    </div>
                                </div>

                                {/* Indications */}
                                <div>
                                    <p className="text-sm text-slate-500 mb-2">Approved Indications</p>
                                    <div className="flex flex-wrap gap-2">
                                        {selectedItem.indications.map((ind, i) => (
                                            <Badge key={i} variant="outline" className="border-slate-700 text-slate-300">
                                                {ind}
                                            </Badge>
                                        ))}
                                    </div>
                                </div>

                                {/* External Links */}
                                <div className="flex gap-3">
                                    <Button variant="outline" size="sm" className="border-slate-700 text-slate-300 hover:bg-slate-800 hover:text-white">
                                        <ExternalLink className="w-4 h-4 mr-2" />
                                        FDA Label
                                    </Button>
                                    <Button variant="outline" size="sm" className="border-slate-700 text-slate-300 hover:bg-slate-800 hover:text-white">
                                        <ExternalLink className="w-4 h-4 mr-2" />
                                        ClinicalTrials.gov
                                    </Button>
                                    <Button variant="outline" size="sm" className="border-slate-700 text-slate-300 hover:bg-slate-800 hover:text-white">
                                        <ExternalLink className="w-4 h-4 mr-2" />
                                        PubMed
                                    </Button>
                                </div>
                            </div>

                            <div className="flex justify-end gap-3">
                                <Button variant="outline" onClick={() => setDetailOpen(false)} className="border-slate-700 text-slate-300 hover:bg-slate-800 hover:text-white">
                                    Close
                                </Button>
                                <Button
                                    className="bg-blue-600 hover:bg-blue-700 text-white"
                                    onClick={() => handleUseAsTemplate(selectedItem)}
                                >
                                    Use this ADC as Template
                                </Button>
                            </div>
                        </>
                    )}
                </DialogContent>
            </Dialog>
        </div>
    )
}
