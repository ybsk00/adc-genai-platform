import { StagingAreaTab } from '@/components/admin/data/StagingAreaTab'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { FlaskConical } from 'lucide-react'

export function StagingAreaPage() {
    return (
        <div className="space-y-6">
            <div>
                <h1 className="text-3xl font-bold text-white flex items-center gap-3">
                    <FlaskConical className="w-8 h-8 text-pink-500" />
                    Staging Area
                </h1>
                <p className="text-slate-400 mt-2">
                    Review AI-refined drafts and promote them to the Golden Set Library.
                </p>
            </div>
            
            <Card className="bg-slate-900 border-slate-800">
                <CardHeader>
                    <CardTitle>Draft Review & Promotion</CardTitle>
                </CardHeader>
                <CardContent>
                    <StagingAreaTab />
                </CardContent>
            </Card>
        </div>
    )
}
