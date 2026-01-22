from apscheduler.schedulers.asyncio import AsyncIOScheduler
from apscheduler.triggers.interval import IntervalTrigger
from apscheduler.triggers.cron import CronTrigger
import logging
from typing import Optional, Callable

logger = logging.getLogger(__name__)

class SchedulerEngine:
    def __init__(self):
        # Timezone 설정 (KST)
        from pytz import timezone
        kst = timezone('Asia/Seoul')
        self.scheduler = AsyncIOScheduler(timezone=kst)
        self.is_started = False

    def start(self):
        """앱 시작 시 스케줄러 가동"""
        if not self.scheduler.running:
            self.scheduler.start()
            self.is_started = True
            
            # 자동 스케줄 등록
            self._register_default_jobs()
            
            logger.info("🚀 Global AsyncIO Scheduler Started! (Timezone: Asia/Seoul)")
            self._log_scheduled_jobs()

    def _log_scheduled_jobs(self):
        """등록된 작업들의 다음 실행 시간 로깅"""
        jobs = self.scheduler.get_jobs()
        for job in jobs:
            logger.info(f"📅 Job [{job.id}] next run: {job.next_run_time}")

    def _register_default_jobs(self):
        """기본 예약 작업 등록"""
        
        # 1. 매일 새벽 4시 Bulk Import (ClinicalTrials.gov 덤프 갱신 반영)
        self.scheduler.add_job(
            self._daily_bulk_import,
            CronTrigger(hour=4, minute=0),
            id="daily_bulk_import",
            replace_existing=True,
            misfire_grace_time=3600
        )
        
        # 2. 매 30분마다 AI Refiner 실행 (소량 배치)
        self.scheduler.add_job(
            self._periodic_refiner,
            IntervalTrigger(minutes=30),
            id="periodic_refiner",
            replace_existing=True,
            misfire_grace_time=1800
        )

    async def _daily_bulk_import(self):
        """새벽 4시 자동 Bulk Import"""
        logger.info("🌅 [Cron] Starting daily bulk import...")
        try:
            from app.services.bulk_importer import BulkImporter
            from app.services.job_lock import job_lock
            
            if not await job_lock.acquire("bulk_import"):
                logger.info("⏭️ [Cron] Bulk import skipped - already running")
                return
            
            try:
                importer = BulkImporter()
                await importer.run_import(max_studies=5000, mode="daily")
            finally:
                await job_lock.release("bulk_import")
                
        except Exception as e:
            logger.error(f"[Cron] Bulk import error: {e}")

    async def _periodic_refiner(self):
        """주기적 AI Refiner (소량 배치)"""
        logger.info("🔄 [Cron] Starting periodic AI refiner...")
        try:
            from app.services.ai_refiner import ai_refiner
            from app.services.job_lock import job_lock
            
            if not await job_lock.acquire("ai_refiner"):
                logger.info("⏭️ [Cron] AI Refiner skipped - already running")
                return
            
            try:
                await ai_refiner.process_pending_records(max_records=20)
            finally:
                await job_lock.release("ai_refiner")
                
        except Exception as e:
            logger.error(f"[Cron] AI Refiner error: {e}")

    def shutdown(self):
        """앱 종료 시 스케줄러 중지"""
        if self.scheduler.running:
            self.scheduler.shutdown()
            logger.info("🛑 Global Scheduler Shutdown.")

    def add_or_update_job(self, job_id: str, func: Callable, hours: int, args: list = None):
        """
        주기적 작업을 등록하거나 업데이트합니다.
        """
        if self.scheduler.get_job(job_id):
            self.scheduler.remove_job(job_id)
            logger.info(f"🔄 Existing job {job_id} removed for update.")

        self.scheduler.add_job(
            func,
            trigger=IntervalTrigger(hours=hours),
            id=job_id,
            args=args or [],
            replace_existing=True,
            misfire_grace_time=3600 # 1시간까지는 지연 실행 허용
        )
        logger.info(f"⏰ Job {job_id} scheduled every {hours} hours.")

    def remove_job(self, job_id: str):
        """작업 제거"""
        if self.scheduler.get_job(job_id):
            self.scheduler.remove_job(job_id)
            logger.info(f"🗑️ Job {job_id} removed.")
    
    def get_scheduled_jobs(self) -> list:
        """등록된 작업 목록"""
        return [
            {"id": job.id, "next_run": str(job.next_run_time)}
            for job in self.scheduler.get_jobs()
        ]

# 전역 인스턴스
scheduler_engine = SchedulerEngine()

