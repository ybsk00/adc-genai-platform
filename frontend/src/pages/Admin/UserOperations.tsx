import { motion } from 'framer-motion'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import {
    Users,
    Activity,
} from 'lucide-react'
import { SimulationLogsTab } from '@/components/admin/users/SimulationLogsTab'
import { UserManagementTab } from '@/components/admin/users/UserManagementTab'
import { Button } from '@/components/ui/button'
import { UserPlus } from 'lucide-react'

export function UserOperations() {
    return (
        <div className="space-y-6">
            {/* Page Header */}
            <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                className="flex items-center justify-between"
            >
                <div>
                    <h1 className="text-2xl font-bold text-white">User Operations</h1>
                    <p className="text-slate-400 mt-1">사용자 관리 및 시뮬레이션 로그 모니터링</p>
                </div>
                <Button className="bg-purple-500 hover:bg-purple-600">
                    <UserPlus className="w-4 h-4 mr-2" />
                    Add User
                </Button>
            </motion.div>

            <Tabs defaultValue="users" className="space-y-4">
                <TabsList className="bg-slate-900 border border-slate-800">
                    <TabsTrigger value="users" className="data-[state=active]:bg-slate-800 text-slate-400 data-[state=active]:text-white">
                        <Users className="w-4 h-4 mr-2" />
                        User Management
                    </TabsTrigger>
                    <TabsTrigger value="logs" className="data-[state=active]:bg-slate-800 text-slate-400 data-[state=active]:text-white">
                        <Activity className="w-4 h-4 mr-2" />
                        Simulation Logs
                    </TabsTrigger>
                </TabsList>

                <TabsContent value="users">
                    <UserManagementTab />
                </TabsContent>

                <TabsContent value="logs">
                    <SimulationLogsTab />
                </TabsContent>
            </Tabs>
        </div>
    )
}
