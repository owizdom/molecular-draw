#!/usr/bin/env python
"""
Celery worker entry point.
Run with: celery -A app.celery_app worker --loglevel=info
Or use this file: celery -A celery_worker worker --loglevel=info
"""
from app.celery_app import celery_app

# Import tasks to register them
from app.tasks import molecular_tasks  # noqa

if __name__ == "__main__":
    celery_app.start()
