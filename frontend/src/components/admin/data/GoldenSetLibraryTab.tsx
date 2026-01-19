import { useState, useEffect, useRef } from 'react'
import {
    useReactTable,
    getCoreRowModel,
    getFilteredRowModel,
    getPaginationRowModel,
    getSortedRowModel,
    flexRender,
} from '@tanstack/react-table'
import type {
    ColumnDef,
    SortingState,
    ColumnFiltersState,
} from '@tanstack/react-table'
import {
    Table,
    TableBody,
    TableCell,
    TableHead,
    TableHeader,
    TableRow,
} from '@/components/ui/table'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import {
    DropdownMenu,
    DropdownMenuCheckboxItem,
    DropdownMenuContent,
    DropdownMenuTrigger,
} from '@/components/ui/dropdown-menu'
import {
    Dialog,
    DialogContent,
    DialogDescription,
    DialogFooter,
    DialogHeader,
    DialogTitle,
    DialogTrigger,
} from "@/components/ui/dialog"
import { Tabs, TabsList, TabsTrigger } from '@/components/ui/tabs'
import { Label } from "@/components/ui/label"
import { Textarea } from '@/components/ui/textarea'
import {
    ChevronDown,
    Check,
    X,
    Loader2,
    RefreshCw,
    Plus,
    Upload,
    FileText,
    AlertCircle,
    Search
} from 'lucide-react'
import { Badge } from '@/components/ui/badge'
import { toast } from 'sonner'
import { format } from 'date-fns'
import SmilesDrawer from 'smiles-drawer'

// API Response Type
export type GoldenSetItem = {
    id: string
    drug_name: string
    target: string
    payload: string | null
    linker: string | null
    enrichment_source: string
    status: 'approved' | 'rejected' | 'draft' | 'ongoing'
    outcome_type?: 'Success' | 'Failure' | 'Terminated'
    failure_reason?: string
    ip_status?: 'Green' | 'Yellow' | 'Red' | 'Unknown'
    smiles_code?: string
    created_at: string
}

// Structure Viewer Component
function StructureViewer({ smiles, id }: { smiles: string, id: string }) {
    const canvasRef = useRef<HTMLCanvasElement>(null)

    useEffect(() => {
        if (canvasRef.current && smiles) {
            const drawer = new SmilesDrawer.Drawer({
                width: 200,
                height: 100,
                compactDrawing: false,
            })
            SmilesDrawer.parse(smiles, (tree: any) => {
                drawer.draw(tree, canvasRef.current, 'light', false)
            }, (err: any) => {
                console.error(err)
            })
        }
    }, [smiles, id])

    return <canvas ref={canvasRef} id={`structure-${id}`} width={200} height={100} />
}

export function GoldenSetLibraryTab() {
    const [data, setData] = useState<GoldenSetItem[]>([])
    const [loading, setLoading] = useState(true)
    const [sorting, setSorting] = useState<SortingState>([])
    const [columnFilters, setColumnFilters] = useState<ColumnFiltersState>([])
    const [rowSelection, setRowSelection] = useState({})
    const [statusFilter, setStatusFilter] = useState<string>('all')

    // Manual Entry State
    const [isAddDialogOpen, setIsAddDialogOpen] = useState(false)
    const [newEntry, setNewEntry] = useState({
        name: '',
        target: '',
        smiles: '',
        description: '',
        outcome: 'Success',
        ipStatus: 'Green'
    })

    const fetchData = async () => {
        setLoading(true)
        try {
            // Mock Data for now as API might not return all fields yet
            // const response = await fetch('/api/library/goldenset')
            // const result = await response.json()

            // Mocking data to match the new requirements
            const mockData: GoldenSetItem[] = [
                {
                    id: '1',
                    drug_name: 'Enhertu',
                    target: 'HER2',
                    payload: 'DXd',
                    linker: 'GGFG',
                    enrichment_source: 'ClinicalTrials.gov',
                    status: 'approved',
                    outcome_type: 'Success',
                    ip_status: 'Green',
                    smiles_code: 'CC1=C(C=C(C=C1)NC(=O)C2(CC2)C(=O)NC3=CC=C(C=C3)C(=O)NC(C)C(=O)NC(CC4=CC=CC=C4)C(=O)N(C)CC(=O)N(C)C(C(C)C)C(=O)O)C', // Trastuzumab deruxtecan payload (approx)
                    created_at: '2025-12-01T10:00:00Z'
                },
                {
                    id: '2',
                    drug_name: 'Trodelvy',
                    target: 'TROP2',
                    payload: 'SN-38',
                    linker: 'CL2A',
                    enrichment_source: 'PubMed',
                    status: 'approved',
                    outcome_type: 'Success',
                    ip_status: 'Green',
                    smiles_code: 'CCC1(C2=C(COC1=O)C(=O)N3CC4=CC5=C(C=C4C3=C2)OCO5)O', // SN-38
                    created_at: '2026-01-10T14:30:00Z'
                },
                {
                    id: '3',
                    drug_name: 'Rovalpituzumab Tesirine',
                    target: 'DLL3',
                    payload: 'SG3249',
                    linker: 'Val-Ala',
                    enrichment_source: 'Manual',
                    status: 'approved', // Actually failed in trials, but approved in library as a record
                    outcome_type: 'Failure',
                    failure_reason: 'Lack of efficacy in Phase 3',
                    ip_status: 'Red',
                    smiles_code: 'CN(C)C(=O)OC1=CC=C(C=C1)C2=CC3=C(C=C2)N(C4=CC=CC=C43)C(=O)C(C)(C)S', // Dummy
                    created_at: '2026-01-15T09:00:00Z'
                }
            ]
            setData(mockData)
        } catch (error) {
            console.error(error)
            toast.error('Failed to load golden set library')
        } finally {
            setLoading(false)
        }
    }

    useEffect(() => {
        fetchData()
    }, [])

    const handleAddEntry = async () => {
        // TODO: API Call
        toast.success('New Golden Set entry added')
        setIsAddDialogOpen(false)
        fetchData()
    }

    const filteredData = data.filter(item => {
        if (statusFilter === 'all') return true
        if (statusFilter === 'approved') return item.outcome_type === 'Success'
        if (statusFilter === 'failed') return item.outcome_type === 'Failure'
        if (statusFilter === 'ongoing') return item.status === 'ongoing'
        return true
    })

    const columns: ColumnDef<GoldenSetItem>[] = [
        {
            accessorKey: 'drug_name',
            header: 'Drug Name',
            cell: ({ row }) => (
                <div>
                    <div className="font-medium text-white">{row.getValue('drug_name')}</div>
                    <div className="text-xs text-slate-500">{row.original.enrichment_source}</div>
                </div>
            ),
        },
        {
            accessorKey: 'target',
            header: 'Target',
            cell: ({ row }) => <div className="text-slate-300">{row.getValue('target')}</div>,
        },
        {
            accessorKey: 'smiles_code',
            header: 'Structure',
            cell: ({ row }) => (
                <div className="w-[100px] h-[50px] bg-white rounded overflow-hidden flex items-center justify-center">
                    {row.original.smiles_code ? (
                        <StructureViewer smiles={row.original.smiles_code} id={row.original.id} />
                    ) : (
                        <span className="text-xs text-slate-400">No Structure</span>
                    )}
                </div>
            ),
        },
        {
            accessorKey: 'outcome_type',
            header: 'Outcome',
            cell: ({ row }) => {
                const outcome = row.original.outcome_type
                return (
                    <Badge variant="outline" className={`
                        ${outcome === 'Success' ? 'border-green-500/30 text-green-400' :
                            outcome === 'Failure' ? 'border-red-500/30 text-red-400' :
                                'border-slate-600 text-slate-400'}
                    `}>
                        {outcome || 'Unknown'}
                    </Badge>
                )
            },
        },
        {
            accessorKey: 'ip_status',
            header: 'IP Status',
            cell: ({ row }) => {
                const status = row.original.ip_status
                return (
                    <div className="flex items-center gap-2">
                        <div className={`w-3 h-3 rounded-full ${status === 'Green' ? 'bg-green-500 shadow-[0_0_8px_rgba(34,197,94,0.5)]' :
                                status === 'Yellow' ? 'bg-yellow-500 shadow-[0_0_8px_rgba(234,179,8,0.5)]' :
                                    status === 'Red' ? 'bg-red-500 shadow-[0_0_8px_rgba(239,68,68,0.5)]' :
                                        'bg-slate-500'
                            }`} />
                        <span className="text-sm text-slate-300">{status}</span>
                    </div>
                )
            },
        },
        {
            id: 'actions',
            header: 'Actions',
            cell: ({ row }) => {
                return (
                    <Button
                        variant="ghost"
                        size="sm"
                        className="text-slate-400 hover:text-white"
                    >
                        Edit
                    </Button>
                )
            }
        }
    ]

    const table = useReactTable({
        data: filteredData,
        columns,
        getCoreRowModel: getCoreRowModel(),
        getPaginationRowModel: getPaginationRowModel(),
        getSortedRowModel: getSortedRowModel(),
        getFilteredRowModel: getFilteredRowModel(),
        onSortingChange: setSorting,
        onColumnFiltersChange: setColumnFilters,
        onRowSelectionChange: setRowSelection,
        state: {
            sorting,
            columnFilters,
            rowSelection,
        },
    })

    return (
        <div className="space-y-4">
            {/* Status Filter Tabs */}
            <div className="flex items-center justify-between">
                <Tabs value={statusFilter} onValueChange={setStatusFilter} className="w-[400px]">
                    <TabsList className="bg-slate-900 border border-slate-800">
                        <TabsTrigger value="all" className="data-[state=active]:bg-slate-800">All</TabsTrigger>
                        <TabsTrigger value="approved" className="data-[state=active]:bg-slate-800 text-green-400">Approved</TabsTrigger>
                        <TabsTrigger value="failed" className="data-[state=active]:bg-slate-800 text-red-400">Failed</TabsTrigger>
                        <TabsTrigger value="ongoing" className="data-[state=active]:bg-slate-800 text-blue-400">Ongoing</TabsTrigger>
                    </TabsList>
                </Tabs>

                <div className="flex gap-2">
                    <Button variant="outline" className="border-slate-700 text-slate-300">
                        <Upload className="w-4 h-4 mr-2" />
                        Bulk Upload
                    </Button>
                    <Dialog open={isAddDialogOpen} onOpenChange={setIsAddDialogOpen}>
                        <DialogTrigger asChild>
                            <Button className="bg-purple-500 hover:bg-purple-600">
                                <Plus className="w-4 h-4 mr-2" />
                                Add Data
                            </Button>
                        </DialogTrigger>
                        <DialogContent className="bg-slate-900 border-slate-800 text-white sm:max-w-[600px]">
                            <DialogHeader>
                                <DialogTitle>Add New Golden Set Entry</DialogTitle>
                                <DialogDescription className="text-slate-400">
                                    Manually add a validated ADC record to the library.
                                </DialogDescription>
                            </DialogHeader>
                            <div className="grid gap-4 py-4">
                                <div className="grid grid-cols-2 gap-4">
                                    <div className="space-y-2">
                                        <Label>Drug Name</Label>
                                        <Input
                                            value={newEntry.name}
                                            onChange={(e) => setNewEntry({ ...newEntry, name: e.target.value })}
                                            className="bg-slate-800 border-slate-700"
                                        />
                                    </div>
                                    <div className="space-y-2">
                                        <Label>Target</Label>
                                        <Input
                                            value={newEntry.target}
                                            onChange={(e) => setNewEntry({ ...newEntry, target: e.target.value })}
                                            className="bg-slate-800 border-slate-700"
                                        />
                                    </div>
                                </div>
                                <div className="space-y-2">
                                    <Label>SMILES Code</Label>
                                    <div className="flex gap-2">
                                        <Input
                                            value={newEntry.smiles}
                                            onChange={(e) => setNewEntry({ ...newEntry, smiles: e.target.value })}
                                            className="bg-slate-800 border-slate-700"
                                            placeholder="e.g. C1=CC=CC=C1"
                                        />
                                    </div>
                                    {/* SMILES Preview */}
                                    <div className="h-[120px] bg-white rounded-lg flex items-center justify-center border border-slate-700 mt-2">
                                        {newEntry.smiles ? (
                                            <StructureViewer smiles={newEntry.smiles} id="preview" />
                                        ) : (
                                            <span className="text-slate-400 text-sm">Structure Preview</span>
                                        )}
                                    </div>
                                </div>
                                <div className="grid grid-cols-2 gap-4">
                                    <div className="space-y-2">
                                        <Label>Outcome</Label>
                                        <select
                                            className="w-full h-10 rounded-md border border-slate-700 bg-slate-800 px-3 py-2 text-sm text-white focus:outline-none focus:ring-2 focus:ring-purple-500"
                                            value={newEntry.outcome}
                                            onChange={(e) => setNewEntry({ ...newEntry, outcome: e.target.value })}
                                        >
                                            <option value="Success">Success</option>
                                            <option value="Failure">Failure</option>
                                            <option value="Terminated">Terminated</option>
                                        </select>
                                    </div>
                                    <div className="space-y-2">
                                        <Label>IP Status</Label>
                                        <select
                                            className="w-full h-10 rounded-md border border-slate-700 bg-slate-800 px-3 py-2 text-sm text-white focus:outline-none focus:ring-2 focus:ring-purple-500"
                                            value={newEntry.ipStatus}
                                            onChange={(e) => setNewEntry({ ...newEntry, ipStatus: e.target.value })}
                                        >
                                            <option value="Green">Green (Safe)</option>
                                            <option value="Yellow">Yellow (Caution)</option>
                                            <option value="Red">Red (Blocked)</option>
                                        </select>
                                    </div>
                                </div>
                            </div>
                            <DialogFooter>
                                <Button variant="outline" onClick={() => setIsAddDialogOpen(false)} className="border-slate-700 text-slate-300">
                                    Cancel
                                </Button>
                                <Button onClick={handleAddEntry} className="bg-purple-500 hover:bg-purple-600">
                                    Save Entry
                                </Button>
                            </DialogFooter>
                        </DialogContent>
                    </Dialog>
                </div>
            </div>

            {/* Table */}
            <div className="rounded-md border border-slate-800">
                <Table>
                    <TableHeader className="bg-slate-900">
                        {table.getHeaderGroups().map((headerGroup) => (
                            <TableRow key={headerGroup.id} className="border-slate-800 hover:bg-transparent">
                                {headerGroup.headers.map((header) => {
                                    return (
                                        <TableHead key={header.id} className="text-slate-400">
                                            {header.isPlaceholder
                                                ? null
                                                : flexRender(
                                                    header.column.columnDef.header,
                                                    header.getContext()
                                                )}
                                        </TableHead>
                                    )
                                })}
                            </TableRow>
                        ))}
                    </TableHeader>
                    <TableBody>
                        {loading ? (
                            <TableRow>
                                <TableCell colSpan={columns.length} className="h-24 text-center">
                                    <div className="flex justify-center items-center text-slate-500">
                                        <Loader2 className="w-6 h-6 animate-spin mr-2" />
                                        Loading library...
                                    </div>
                                </TableCell>
                            </TableRow>
                        ) : table.getRowModel().rows?.length ? (
                            table.getRowModel().rows.map((row) => (
                                <TableRow
                                    key={row.id}
                                    data-state={row.getIsSelected() && 'selected'}
                                    className="border-slate-800 hover:bg-slate-800/50"
                                >
                                    {row.getVisibleCells().map((cell) => (
                                        <TableCell key={cell.id}>
                                            {flexRender(
                                                cell.column.columnDef.cell,
                                                cell.getContext()
                                            )}
                                        </TableCell>
                                    ))}
                                </TableRow>
                            ))
                        ) : (
                            <TableRow>
                                <TableCell
                                    colSpan={columns.length}
                                    className="h-24 text-center text-slate-500"
                                >
                                    No records found.
                                </TableCell>
                            </TableRow>
                        )}
                    </TableBody>
                </Table>
            </div>
        </div>
    )
}
