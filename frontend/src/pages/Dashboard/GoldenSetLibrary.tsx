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
                <h1 className="text-2xl font-bold text-gray-900">Golden Set Library</h1>
                <p className="text-gray-500 mt-1">FDA 승인 ADC 데이터베이스 - 15개의 검증된 ADC 정보</p>
            </motion.div>

            {/* Search & Filters */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.1 }}
            >
                <Card>
                    <CardContent className="p-4">
                        <div className="flex flex-col md:flex-row gap-4">
                            {/* Search */}
                            <div className="relative flex-1">
                                <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-gray-400" />
                                <Input
                                    placeholder="약물명, 회사명으로 검색..."
                                    value={searchQuery}
                                    onChange={(e) => setSearchQuery(e.target.value)}
                                    className="pl-10"
                                />
                            </div>

                            {/* Target Filter */}
                            <Select value={targetFilter} onValueChange={setTargetFilter}>
                                <SelectTrigger className="w-full md:w-[180px]">
                                    <Target className="w-4 h-4 mr-2" />
                                    <SelectValue placeholder="타겟" />
                                </SelectTrigger>
                                <SelectContent>
                                    <SelectItem value="all">모든 타겟</SelectItem>
                                    {targets.map(t => (
                                        <SelectItem key={t} value={t}>{t}</SelectItem>
                                    ))}
                                </SelectContent>
                            </Select>

                            {/* Payload Filter */}
                            <Select value={payloadFilter} onValueChange={setPayloadFilter}>
                                <SelectTrigger className="w-full md:w-[180px]">
                                    <Pill className="w-4 h-4 mr-2" />
                                    <SelectValue placeholder="페이로드" />
                                </SelectTrigger>
                                <SelectContent>
                                    <SelectItem value="all">모든 페이로드</SelectItem>
                                    {payloads.map(p => (
                                        <SelectItem key={p} value={p}>{p}</SelectItem>
                                    ))}
                                </SelectContent>
                            </Select>
                        </div>
                    </CardContent>
                </Card>
            </motion.div>

            {/* Results Count */}
            <div className="flex items-center justify-between">
                <p className="text-sm text-gray-500">
                    {filteredData.length}개의 결과
                </p>
                <Button variant="outline" size="sm">
                    <Filter className="w-4 h-4 mr-2" />
                    고급 필터
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
                        <Card className="hover:shadow-lg transition-shadow cursor-pointer group">
                            <CardHeader className="pb-3">
                                <div className="flex items-start justify-between">
                                    <div>
                                        <CardTitle className="text-lg group-hover:text-[#007AFF] transition-colors">
                                            {item.drugName}
                                        </CardTitle>
                                        <CardDescription className="text-xs mt-1">
                                            {item.brandName}
                                        </CardDescription>
                                    </div>
                                    <Badge className="bg-green-100 text-green-700 hover:bg-green-100">
                                        {item.status}
                                    </Badge>
                                </div>
                            </CardHeader>
                            <CardContent className="space-y-3">
                                <div className="grid grid-cols-2 gap-2 text-sm">
                                    <div className="flex items-center gap-2 text-gray-600">
                                        <Target className="w-4 h-4 text-[#007AFF]" />
                                        <span>{item.target}</span>
                                    </div>
                                    <div className="flex items-center gap-2 text-gray-600">
                                        <Pill className="w-4 h-4 text-purple-500" />
                                        <span>{item.payload}</span>
                                    </div>
                                    <div className="flex items-center gap-2 text-gray-600">
                                        <Link2 className="w-4 h-4 text-orange-500" />
                                        <span className="truncate">{item.linker.split(' ')[0]}</span>
                                    </div>
                                    <div className="flex items-center gap-2 text-gray-600">
                                        <FlaskConical className="w-4 h-4 text-green-500" />
                                        <span>DAR {item.dar}</span>
                                    </div>
                                </div>

                                <div className="flex items-center gap-2 text-sm text-gray-500">
                                    <Building2 className="w-4 h-4" />
                                    <span className="truncate">{item.company}</span>
                                </div>

                                <div className="flex items-center gap-2 text-sm text-gray-500">
                                    <Calendar className="w-4 h-4" />
                                    <span>승인일: {item.approvalDate}</span>
                                </div>

                                <div className="flex gap-2 pt-2">
                                    <Button
                                        variant="outline"
                                        size="sm"
                                        className="flex-1"
                                        onClick={() => handleViewDetail(item)}
                                    >
                                        상세보기
                                    </Button>
                                    <Button
                                        size="sm"
                                        className="flex-1 bg-[#007AFF] hover:bg-[#0066DD]"
                                        onClick={() => handleUseAsTemplate(item)}
                                    >
                                        템플릿 사용
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
                <DialogContent className="max-w-2xl">
                    {selectedItem && (
                        <>
                            <DialogHeader>
                                <DialogTitle className="text-xl">{selectedItem.drugName}</DialogTitle>
                                <DialogDescription>{selectedItem.brandName}</DialogDescription>
                            </DialogHeader>

                            <div className="space-y-6 py-4">
                                {/* Basic Info */}
                                <div className="grid grid-cols-2 gap-4">
                                    <div>
                                        <p className="text-sm text-gray-500">타겟</p>
                                        <p className="font-medium">{selectedItem.target}</p>
                                    </div>
                                    <div>
                                        <p className="text-sm text-gray-500">페이로드</p>
                                        <p className="font-medium">{selectedItem.payload}</p>
                                    </div>
                                    <div>
                                        <p className="text-sm text-gray-500">항체</p>
                                        <p className="font-medium">{selectedItem.antibody}</p>
                                    </div>
                                    <div>
                                        <p className="text-sm text-gray-500">링커</p>
                                        <p className="font-medium">{selectedItem.linker}</p>
                                    </div>
                                    <div>
                                        <p className="text-sm text-gray-500">DAR</p>
                                        <p className="font-medium">{selectedItem.dar}</p>
                                    </div>
                                    <div>
                                        <p className="text-sm text-gray-500">개발사</p>
                                        <p className="font-medium">{selectedItem.company}</p>
                                    </div>
                                </div>

                                {/* Indications */}
                                <div>
                                    <p className="text-sm text-gray-500 mb-2">승인 적응증</p>
                                    <div className="flex flex-wrap gap-2">
                                        {selectedItem.indications.map((ind, i) => (
                                            <Badge key={i} variant="outline">{ind}</Badge>
                                        ))}
                                    </div>
                                </div>

                                {/* External Links */}
                                <div className="flex gap-3">
                                    <Button variant="outline" size="sm">
                                        <ExternalLink className="w-4 h-4 mr-2" />
                                        FDA Label
                                    </Button>
                                    <Button variant="outline" size="sm">
                                        <ExternalLink className="w-4 h-4 mr-2" />
                                        ClinicalTrials.gov
                                    </Button>
                                    <Button variant="outline" size="sm">
                                        <ExternalLink className="w-4 h-4 mr-2" />
                                        PubMed
                                    </Button>
                                </div>
                            </div>

                            <div className="flex justify-end gap-3">
                                <Button variant="outline" onClick={() => setDetailOpen(false)}>
                                    닫기
                                </Button>
                                <Button
                                    className="bg-[#007AFF] hover:bg-[#0066DD]"
                                    onClick={() => handleUseAsTemplate(selectedItem)}
                                >
                                    이 ADC를 템플릿으로 사용
                                </Button>
                            </div>
                        </>
                    )}
                </DialogContent>
            </Dialog>
        </div>
    )
}
