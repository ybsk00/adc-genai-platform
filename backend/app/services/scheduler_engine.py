from apscheduler.schedulers.asyncio import AsyncIOScheduler
from apscheduler.triggers.interval import IntervalTrigger
import logging
from typing import Optional, Callable

logger = logging.getLogger(__name__)

class SchedulerEngine:
    def __init__(self):
        self.scheduler = AsyncIOScheduler()
        self.is_started = False

    def start(self):
        """앱 시작 시 스케줄러 가동"""
        if not self.scheduler.running:
            self.scheduler.start()
            self.is_started = True
            logger.info("🚀 Global AsyncIO Scheduler Started!")

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

# 전역 인스턴스
scheduler_engine = SchedulerEngine()
