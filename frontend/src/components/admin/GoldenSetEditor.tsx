import { useState } from 'react'
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
import { ChevronDown, Plus, Upload, Save, Trash2 } from 'lucide-react'
import { Badge } from '@/components/ui/badge'
import { toast } from 'sonner'

// Mock Data Type
export type GoldenSet = {
    id: string
    drugName: string
    target: string
    antibody: string
    payload: string
    phase: 'Approved' | 'Phase 3' | 'Phase 2' | 'Phase 1' | 'Pre-clinical' | 'Withdrawn'
    status: 'Active' | 'Discontinued' | 'Withdrawn'
}

// Initial Mock Data
const initialData: GoldenSet[] = [
    { id: '1', drugName: 'Enhertu', target: 'HER2', antibody: 'Trastuzumab', payload: 'DXd', phase: 'Approved', status: 'Active' },
    { id: '2', drugName: 'Kadcyla', target: 'HER2', antibody: 'Trastuzumab', payload: 'DM1', phase: 'Approved', status: 'Active' },
    { id: '3', drugName: 'Padcev', target: 'Nectin-4', antibody: 'Enfortumab', payload: 'MMAE', phase: 'Approved', status: 'Active' },
    { id: '4', drugName: 'Trodelvy', target: 'TROP-2', antibody: 'Sacituzumab', payload: 'SN-38', phase: 'Approved', status: 'Active' },
    { id: '5', drugName: 'Blenrep', target: 'BCMA', antibody: 'Belantamab', payload: 'MMAF', phase: 'Withdrawn', status: 'Discontinued' },
]

export function GoldenSetEditor() {
    const [data, setData] = useState<GoldenSet[]>(initialData)
    const [sorting, setSorting] = useState<SortingState>([])
    const [columnFilters, setColumnFilters] = useState<ColumnFiltersState>([])
    const [rowSelection, setRowSelection] = useState({})

    // Editable Cell Component
    const EditableCell = ({ getValue, row, column, table }: any) => {
        const initialValue = getValue()
        const [value, setValue] = useState(initialValue)

        const onBlur = () => {
            table.options.meta?.updateData(row.index, column.id, value)
        }

        return (
            <Input
                value={value}
                onChange={e => setValue(e.target.value)}
                onBlur={onBlur}
                className="h-8 border-transparent hover:border-slate-700 bg-transparent focus:bg-slate-800 text-slate-300"
            />
        )
    }

    const columns: ColumnDef<GoldenSet>[] = [
        {
            accessorKey: 'drugName',
            header: 'Drug Name',
            cell: EditableCell,
        },
        {
            accessorKey: 'target',
            header: 'Target',
            cell: EditableCell,
        },
        {
            accessorKey: 'antibody',
            header: 'Antibody',
            cell: EditableCell,
        },
        {
            accessorKey: 'payload',
            header: 'Payload',
            cell: EditableCell,
        },
        {
            accessorKey: 'phase',
            header: 'Phase',
            cell: ({ row }) => {
                const phase = row.getValue('phase') as string
                return (
                    <Badge variant="outline" className={
                        phase === 'Approved' ? 'border-green-500/30 text-green-400' :
                            phase.includes('Phase') ? 'border-blue-500/30 text-blue-400' :
                                'border-slate-600 text-slate-400'
                    }>
                        {phase}
                    </Badge>
                )
            }
        },
        {
            accessorKey: 'status',
            header: 'Status',
            cell: ({ row }) => {
                const status = row.getValue('status') as string
                return (
                    <span className={status === 'Active' ? 'text-green-400' : 'text-red-400'}>
                        {status}
                    </span>
                )
            }
        },
        {
            id: 'actions',
            cell: ({ row }) => {
                return (
                    <Button
                        variant="ghost"
                        size="icon"
                        onClick={() => {
                            const newData = [...data]
                            newData.splice(row.index, 1)
                            setData(newData)
                            toast.success('Row deleted')
                        }}
                        className="text-slate-500 hover:text-red-400"
                    >
                        <Trash2 className="w-4 h-4" />
                    </Button>
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
        meta: {
            updateData: (rowIndex: number, columnId: string, value: unknown) => {
                setData(old =>
                    old.map((row, index) => {
                        if (index === rowIndex) {
                            return {
                                ...old[rowIndex]!,
                                [columnId]: value,
                            }
                        }
                        return row
                    })
                )
            },
        },
    })

    const handleAddRow = () => {
        const newRow: GoldenSet = {
            id: Math.random().toString(36).substr(2, 9),
            drugName: 'New Drug',
            target: '-',
            antibody: '-',
            payload: '-',
            phase: 'Pre-clinical',
            status: 'Active'
        }
        setData([...data, newRow])
        toast.success('New row added')
    }

    const handleSave = () => {
        // TODO: API Call to save data
        toast.success('Changes saved to database')
    }

    return (
        <div className="space-y-4">
            <div className="flex items-center justify-between">
                <div className="flex items-center gap-2">
                    <Input
                        placeholder="Filter drugs..."
                        value={(table.getColumn('drugName')?.getFilterValue() as string) ?? ''}
                        onChange={(event) =>
                            table.getColumn('drugName')?.setFilterValue(event.target.value)
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
                                            {column.id}
                                        </DropdownMenuCheckboxItem>
                                    )
                                })}
                        </DropdownMenuContent>
                    </DropdownMenu>
                </div>
                <div className="flex gap-2">
                    <Button variant="outline" className="border-slate-700 text-slate-300" onClick={() => toast.info('CSV Import feature coming soon')}>
                        <Upload className="w-4 h-4 mr-2" />
                        Import CSV
                    </Button>
                    <Button onClick={handleAddRow} className="bg-slate-700 hover:bg-slate-600">
                        <Plus className="w-4 h-4 mr-2" />
                        Add Row
                    </Button>
                    <Button onClick={handleSave} className="bg-purple-600 hover:bg-purple-700">
                        <Save className="w-4 h-4 mr-2" />
                        Save Changes
                    </Button>
                </div>
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
                        {table.getRowModel().rows?.length ? (
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
                                    No results.
                                </TableCell>
                            </TableRow>
                        )}
                    </TableBody>
                </Table>
            </div>
            <div className="flex items-center justify-end space-x-2 py-4">
                <div className="flex-1 text-sm text-slate-500">
                    {table.getFilteredSelectedRowModel().rows.length} of{' '}
                    {table.getFilteredRowModel().rows.length} row(s) selected.
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
        </div>
    )
}
