/**
 * BOM Table - Bill of Materials 테이블
 * commercial_reagents DB 연동 시약 리스트
 */
import { useState } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow
} from '@/components/ui/table'
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger
} from '@/components/ui/tooltip'
import { ScrollArea } from '@/components/ui/scroll-area'
import { motion } from 'framer-motion'
import {
  Package,
  DollarSign,
  ExternalLink,
  Search,
  AlertTriangle,
  CheckCircle2,
  Clock,
  Building2,
  FlaskConical,
  ArrowUpDown
} from 'lucide-react'

export interface ReagentItem {
  id: string
  name: string
  casNumber?: string
  smiles?: string
  supplier: string
  supplierUrl?: string
  catalogNumber?: string
  quantity: string
  unitPrice: number
  currency: string
  availability: 'in-stock' | 'limited' | 'out-of-stock' | 'custom-synthesis'
  leadTime?: string
  purity?: string
  notes?: string
}

interface BOMTableProps {
  items: ReagentItem[]
  title?: string
  onSupplierClick?: (item: ReagentItem) => void
  className?: string
}

export function BOMTable({
  items,
  title = 'Bill of Materials',
  onSupplierClick,
  className
}: BOMTableProps) {
  const [searchTerm, setSearchTerm] = useState('')
  const [sortBy, setSortBy] = useState<'price' | 'name' | 'supplier'>('name')
  const [sortOrder, setSortOrder] = useState<'asc' | 'desc'>('asc')

  // Filter and sort items
  const filteredItems = items
    .filter(item =>
      item.name.toLowerCase().includes(searchTerm.toLowerCase()) ||
      item.supplier.toLowerCase().includes(searchTerm.toLowerCase()) ||
      item.casNumber?.includes(searchTerm)
    )
    .sort((a, b) => {
      let comparison = 0
      switch (sortBy) {
        case 'price':
          comparison = a.unitPrice - b.unitPrice
          break
        case 'name':
          comparison = a.name.localeCompare(b.name)
          break
        case 'supplier':
          comparison = a.supplier.localeCompare(b.supplier)
          break
      }
      return sortOrder === 'asc' ? comparison : -comparison
    })

  // Calculate totals
  const totalCost = items.reduce((sum, item) => sum + item.unitPrice, 0)
  const inStockCount = items.filter(i => i.availability === 'in-stock').length
  const limitedCount = items.filter(i => i.availability === 'limited').length
  const outOfStockCount = items.filter(i =>
    i.availability === 'out-of-stock' || i.availability === 'custom-synthesis'
  ).length

  const handleSort = (field: typeof sortBy) => {
    if (sortBy === field) {
      setSortOrder(sortOrder === 'asc' ? 'desc' : 'asc')
    } else {
      setSortBy(field)
      setSortOrder('asc')
    }
  }

  const getAvailabilityConfig = (availability: ReagentItem['availability']) => {
    switch (availability) {
      case 'in-stock':
        return {
          icon: <CheckCircle2 className="w-3 h-3" />,
          color: 'text-green-600 bg-green-50',
          label: 'In Stock'
        }
      case 'limited':
        return {
          icon: <AlertTriangle className="w-3 h-3" />,
          color: 'text-amber-600 bg-amber-50',
          label: 'Limited'
        }
      case 'out-of-stock':
        return {
          icon: <Clock className="w-3 h-3" />,
          color: 'text-red-600 bg-red-50',
          label: 'Out of Stock'
        }
      case 'custom-synthesis':
        return {
          icon: <FlaskConical className="w-3 h-3" />,
          color: 'text-purple-600 bg-purple-50',
          label: 'Custom'
        }
    }
  }

  return (
    <Card className={className}>
      <CardHeader className="pb-2">
        <div className="flex items-center justify-between">
          <CardTitle className="text-sm font-medium flex items-center gap-2">
            <Package className="w-4 h-4 text-orange-500" />
            {title}
          </CardTitle>
          <Badge variant="outline" className="text-xs">
            {items.length} items
          </Badge>
        </div>

        {/* Summary Stats */}
        <div className="flex flex-wrap gap-2 mt-2">
          <Badge variant="secondary" className="text-[10px]">
            <DollarSign className="w-3 h-3 mr-1" />
            Total: ${totalCost.toLocaleString()}
          </Badge>
          <Badge variant="secondary" className="text-[10px] bg-green-50 text-green-600">
            <CheckCircle2 className="w-3 h-3 mr-1" />
            {inStockCount} In Stock
          </Badge>
          {limitedCount > 0 && (
            <Badge variant="secondary" className="text-[10px] bg-amber-50 text-amber-600">
              <AlertTriangle className="w-3 h-3 mr-1" />
              {limitedCount} Limited
            </Badge>
          )}
          {outOfStockCount > 0 && (
            <Badge variant="secondary" className="text-[10px] bg-red-50 text-red-600">
              <Clock className="w-3 h-3 mr-1" />
              {outOfStockCount} Unavailable
            </Badge>
          )}
        </div>

        {/* Search */}
        <div className="relative mt-2">
          <Search className="absolute left-2 top-1/2 -translate-y-1/2 w-4 h-4 text-gray-400" />
          <Input
            placeholder="Search reagents, CAS numbers, suppliers..."
            value={searchTerm}
            onChange={(e) => setSearchTerm(e.target.value)}
            className="pl-8 h-8 text-xs"
          />
        </div>
      </CardHeader>

      <CardContent className="p-0">
        <ScrollArea className="h-[350px]">
          <Table>
            <TableHeader className="sticky top-0 bg-white z-10">
              <TableRow>
                <TableHead className="w-[200px]">
                  <Button
                    variant="ghost"
                    size="sm"
                    className="h-6 text-xs font-medium"
                    onClick={() => handleSort('name')}
                  >
                    Reagent
                    <ArrowUpDown className="w-3 h-3 ml-1" />
                  </Button>
                </TableHead>
                <TableHead>
                  <Button
                    variant="ghost"
                    size="sm"
                    className="h-6 text-xs font-medium"
                    onClick={() => handleSort('supplier')}
                  >
                    Supplier
                    <ArrowUpDown className="w-3 h-3 ml-1" />
                  </Button>
                </TableHead>
                <TableHead className="text-right">
                  <Button
                    variant="ghost"
                    size="sm"
                    className="h-6 text-xs font-medium"
                    onClick={() => handleSort('price')}
                  >
                    Price
                    <ArrowUpDown className="w-3 h-3 ml-1" />
                  </Button>
                </TableHead>
                <TableHead className="text-center">Status</TableHead>
              </TableRow>
            </TableHeader>
            <TableBody>
              {filteredItems.map((item, index) => (
                <motion.tr
                  key={item.id}
                  initial={{ opacity: 0, y: 10 }}
                  animate={{ opacity: 1, y: 0 }}
                  transition={{ delay: index * 0.02 }}
                  className="group hover:bg-slate-50"
                >
                  <TableCell className="py-2">
                    <div>
                      <p className="text-xs font-medium truncate max-w-[180px]">
                        {item.name}
                      </p>
                      {item.casNumber && (
                        <p className="text-[10px] text-gray-500 font-mono">
                          CAS: {item.casNumber}
                        </p>
                      )}
                      <p className="text-[10px] text-gray-400">
                        {item.quantity}
                        {item.purity && ` | ${item.purity}`}
                      </p>
                    </div>
                  </TableCell>
                  <TableCell className="py-2">
                    <TooltipProvider>
                      <Tooltip>
                        <TooltipTrigger asChild>
                          <div
                            className="flex items-center gap-1 cursor-pointer hover:text-blue-600"
                            onClick={() => onSupplierClick?.(item)}
                          >
                            <Building2 className="w-3 h-3 text-gray-400" />
                            <span className="text-xs truncate max-w-[100px]">
                              {item.supplier}
                            </span>
                            {item.supplierUrl && (
                              <ExternalLink className="w-3 h-3 text-gray-400 opacity-0 group-hover:opacity-100" />
                            )}
                          </div>
                        </TooltipTrigger>
                        <TooltipContent>
                          <p className="text-xs">
                            {item.catalogNumber && `Cat# ${item.catalogNumber}`}
                            {item.leadTime && ` | Lead time: ${item.leadTime}`}
                          </p>
                        </TooltipContent>
                      </Tooltip>
                    </TooltipProvider>
                  </TableCell>
                  <TableCell className="text-right py-2">
                    <span className="text-xs font-semibold text-gray-800">
                      {item.currency === 'USD' ? '$' : item.currency}
                      {item.unitPrice.toLocaleString()}
                    </span>
                  </TableCell>
                  <TableCell className="text-center py-2">
                    <Badge
                      variant="secondary"
                      className={`text-[10px] ${
                        getAvailabilityConfig(item.availability).color
                      }`}
                    >
                      {getAvailabilityConfig(item.availability).icon}
                      <span className="ml-1">
                        {getAvailabilityConfig(item.availability).label}
                      </span>
                    </Badge>
                  </TableCell>
                </motion.tr>
              ))}
            </TableBody>
          </Table>

          {filteredItems.length === 0 && (
            <div className="text-center py-8 text-gray-400">
              <Package className="w-8 h-8 mx-auto mb-2 opacity-50" />
              <p className="text-sm">No reagents found</p>
            </div>
          )}
        </ScrollArea>
      </CardContent>
    </Card>
  )
}

/**
 * Compact BOM Summary
 */
export function BOMSummary({
  totalItems,
  totalCost,
  availableCount,
  currency = 'USD'
}: {
  totalItems: number
  totalCost: number
  availableCount: number
  currency?: string
}) {
  const availablePercent = (availableCount / totalItems) * 100

  return (
    <div className="flex items-center gap-3 p-3 rounded-lg bg-orange-50 border border-orange-200">
      <Package className="w-5 h-5 text-orange-600" />
      <div className="flex-1">
        <div className="flex items-center justify-between">
          <span className="text-xs text-orange-800">
            {totalItems} reagents | {availableCount} available
          </span>
          <span className="text-sm font-semibold text-orange-700">
            {currency === 'USD' ? '$' : currency}{totalCost.toLocaleString()}
          </span>
        </div>
        <div className="mt-1 h-1.5 bg-orange-100 rounded-full overflow-hidden">
          <div
            className="h-full bg-orange-500 rounded-full"
            style={{ width: `${availablePercent}%` }}
          />
        </div>
      </div>
    </div>
  )
}
