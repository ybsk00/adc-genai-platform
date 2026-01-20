import { useState, useEffect } from 'react'
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
} from "@/components/ui/dialog"
import { Label } from "@/components/ui/label"
import { ChevronDown, Check, X, Loader2, RefreshCw, AlertCircle, Clock, CheckCircle2, Sparkles } from 'lucide-react'
import { Badge } from '@/components/ui/badge'
import {
    Tooltip,
    TooltipContent,
    TooltipProvider,
    TooltipTrigger,
} from "@/components/ui/tooltip"
import { toast } from 'sonner'
import { format } from 'date-fns'

// API Response Type
export type GoldenSetDraft = {
    id: string
    drug_name: string
    target: string
    payload: string | null
    linker: string | null
    enrichment_source: string
    created_at: string
    outcome_type?: string
    failure_reason?: string
}

export function GoldenSetEditor() {
    const [data, setData] = useState<GoldenSetDraft[]>([])
    const [loading, setLoading] = useState(true)
    const [sorting, setSorting] = useState<SortingState>([])
    const [columnFilters, setColumnFilters] = useState<ColumnFiltersState>([])
    const [rowSelection, setRowSelection] = useState({})

    // Reject Dialog State
    const [rejectDialogOpen, setRejectDialogOpen] = useState(false)
    const [selectedId, setSelectedId] = useState<string | null>(null)
    const [rejectReason, setRejectReason] = useState("")

    const fetchData = async () => {
        setLoading(true)
        try {
            const response = await fetch('/api/admin/goldenset/drafts')
            if (!response.ok) throw new Error('Failed to fetch drafts')
            const result = await response.json()
            setData(result)
        } catch (error) {
            console.error(error)
            toast.error('Failed to load golden set drafts')
        } finally {
            setLoading(false)
        }
    }

    useEffect(() => {
        fetchData()
    }, [])

    const handleApprove = async (id: string) => {
        try {
            const response = await fetch(`/api/admin/goldenset/${id}/approve`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({})
            })

            if (!response.ok) throw new Error('Failed to approve')

            toast.success('Golden Set approved and indexed')
            fetchData() // Refresh list
        } catch (error) {
            console.error(error)
            toast.error('Failed to approve golden set')
        }
    }

    const handleReject = async () => {
        if (!selectedId || !rejectReason) return

        try {
            const response = await fetch(`/api/admin/goldenset/${selectedId}/reject`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ reason: rejectReason })
            })

            if (!response.ok) throw new Error('Failed to reject')

            toast.success('Golden Set rejected')
            setRejectDialogOpen(false)
            setRejectReason("")
            setSelectedId(null)
            fetchData() // Refresh list
        } catch (error) {
            console.error(error)
            toast.error('Failed to reject golden set')
        }
    }

    const openRejectDialog = (id: string) => {
        setSelectedId(id)
        setRejectDialogOpen(true)
    }

    const columns: ColumnDef<GoldenSetDraft>[] = [
        accessorKey: 'drug_name',
        header: 'Drug Name',
        cell: ({ row }) => (
            <div className="flex items-center gap-2">
                <span className="font-medium text-white">{row.getValue('drug_name')}</span>
                {row.original.enrichment_source === 'clinical_trials_ai_refined' && (
                    <Badge variant="secondary" className="bg-purple-500/10 text-purple-400 border-purple-500/20 h-5 px-1.5 text-[10px] gap-1">
                        <Sparkles className="w-2.5 h-2.5" /> AI
                    </Badge>
                )}
            </div>
        ),
        },
{
    accessorKey: 'outcome_type',
        header: 'Outcome',
            cell: ({ row }) => {
                const type = row.original.outcome_type
                const reason = row.original.failure_reason

                if (type === 'Success') {
                    return (
                        <Badge variant="outline" className="bg-green-500/10 text-green-400 border-green-500/20 flex w-fit gap-1 items-center">
                            <CheckCircle2 className="w-3 h-3" /> Success
                        </Badge>
                    )
                }

                if (type === 'Failure') {
                    return (
                        <TooltipProvider>
                            <Tooltip>
                                <TooltipTrigger>
                                    <Badge variant="outline" className="bg-red-500/10 text-red-400 border-red-500/20 flex w-fit gap-1 items-center cursor-help">
                                        <AlertCircle className="w-3 h-3" /> Failure
                                    </Badge>
                                </TooltipTrigger>
                                <TooltipContent className="bg-slate-900 border-slate-700 text-slate-300 max-w-xs">
                                    <p>Reason: {reason || "Unknown"}</p>
                                </TooltipContent>
                            </Tooltip>
                        </TooltipProvider>
                    )
                }

                if (type === 'Ongoing') {
                    return (
                        <Badge variant="outline" className="bg-blue-500/10 text-blue-400 border-blue-500/20 flex w-fit gap-1 items-center">
                            <Clock className="w-3 h-3" /> Ongoing
                        </Badge>
                    )
                }

                return <span className="text-slate-500 text-xs">-</span>
            }
},
{
    accessorKey: 'target',
        header: 'Target',
            cell: ({ row }) => <div className="text-slate-300">{row.getValue('target')}</div>,
        },
{
    accessorKey: 'payload',
        header: 'Payload',
            cell: ({ row }) => <div className="text-slate-300">{row.getValue('payload') || '-'}</div>,
        },
{
    accessorKey: 'linker',
        header: 'Linker',
            cell: ({ row }) => <div className="text-slate-300">{row.getValue('linker') || '-'}</div>,
        },
{
    accessorKey: 'enrichment_source',
        header: 'Source',
            cell: ({ row }) => (
                <Badge variant="outline" className="border-slate-700 text-slate-400">
                    {row.getValue('enrichment_source')}
                </Badge>
            ),
        },
{
    accessorKey: 'created_at',
        header: 'Date',
            cell: ({ row }) => {
                const date = new Date(row.getValue('created_at'))
                return <div className="text-slate-400 text-xs">{format(date, 'yyyy-MM-dd HH:mm')}</div>
            },
        },
{
    id: 'actions',
        header: 'Actions',
            cell: ({ row }) => {
                return (
                    <div className="flex items-center gap-2">
                        <Button
                            variant="ghost"
                            size="icon"
                            onClick={() => handleApprove(row.original.id)}
                            className="text-green-500 hover:text-green-400 hover:bg-green-500/10"
                            title="Approve"
                        >
                            <Check className="w-4 h-4" />
                        </Button>
                        <Button
                            variant="ghost"
                            size="icon"
                            onClick={() => openRejectDialog(row.original.id)}
                            className="text-red-500 hover:text-red-400 hover:bg-red-500/10"
                            title="Reject"
                        >
                            <X className="w-4 h-4" />
                        </Button>
                    </div>
                )
            }
}
    ]

const table = useReactTable({
    data,
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
        <div className="flex items-center justify-between">
            <div className="flex items-center gap-2">
                <Input
                    placeholder="Filter drugs..."
                    value={(table.getColumn('drug_name')?.getFilterValue() as string) ?? ''}
                    onChange={(event) =>
                        table.getColumn('drug_name')?.setFilterValue(event.target.value)
                    }
                    className="max-w-sm bg-slate-800 border-slate-700 text-white"
                />
                <DropdownMenu>
                    <DropdownMenuTrigger asChild>
                        <Button variant="outline" className="ml-auto border-slate-700 text-slate-300">
                            Columns <ChevronDown className="ml-2 h-4 w-4" />
                        </Button>
                    </DropdownMenuTrigger>
                    <DropdownMenuContent align="end" className="bg-slate-800 border-slate-700">
                        {table
                            .getAllColumns()
                            .filter((column) => column.getCanHide())
                            .map((column) => {
                                return (
                                    <DropdownMenuCheckboxItem
                                        key={column.id}
                                        className="capitalize text-slate-300 focus:bg-slate-700"
                                        checked={column.getIsVisible()}
                                        onCheckedChange={(value) =>
                                            column.toggleVisibility(!!value)
                                        }
                                    >
                                        {column.id.replace('_', ' ')}
                                    </DropdownMenuCheckboxItem>
                                )
                            })}
                    </DropdownMenuContent>
                </DropdownMenu>
            </div>
            <Button
                variant="outline"
                onClick={fetchData}
                disabled={loading}
                className="border-slate-700 text-slate-300"
            >
                <RefreshCw className={`w-4 h-4 mr-2 ${loading ? 'animate-spin' : ''}`} />
                Refresh
            </Button>
        </div>

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
                                    Loading drafts...
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
                                No pending drafts found.
                            </TableCell>
                        </TableRow>
                    )}
                </TableBody>
            </Table>
        </div>

        <div className="flex items-center justify-end space-x-2 py-4">
            <div className="flex-1 text-sm text-slate-500">
                {table.getFilteredRowModel().rows.length} row(s)
            </div>
            <div className="space-x-2">
                <Button
                    variant="outline"
                    size="sm"
                    onClick={() => table.previousPage()}
                    disabled={!table.getCanPreviousPage()}
                    className="border-slate-700 text-slate-300"
                >
                    Previous
                </Button>
                <Button
                    variant="outline"
                    size="sm"
                    onClick={() => table.nextPage()}
                    disabled={!table.getCanNextPage()}
                    className="border-slate-700 text-slate-300"
                >
                    Next
                </Button>
            </div>
        </div>

        <Dialog open={rejectDialogOpen} onOpenChange={setRejectDialogOpen}>
            <DialogContent className="bg-slate-900 border-slate-800 text-white">
                <DialogHeader>
                    <DialogTitle>Reject Golden Set Draft</DialogTitle>
                    <DialogDescription className="text-slate-400">
                        Please provide a reason for rejecting this data.
                    </DialogDescription>
                </DialogHeader>
                <div className="grid gap-4 py-4">
                    <div className="grid gap-2">
                        <Label htmlFor="reason" className="text-slate-300">Reason</Label>
                        <Input
                            id="reason"
                            value={rejectReason}
                            onChange={(e) => setRejectReason(e.target.value)}
                            className="bg-slate-800 border-slate-700 text-white"
                            placeholder="e.g., Duplicate entry, Incorrect data..."
                        />
                    </div>
                </div>
                <DialogFooter>
                    <Button variant="outline" onClick={() => setRejectDialogOpen(false)} className="border-slate-700 text-slate-300">
                        Cancel
                    </Button>
                    <Button onClick={handleReject} variant="destructive">
                        Reject
                    </Button>
                </DialogFooter>
            </DialogContent>
        </Dialog>
    </div>
)
}
