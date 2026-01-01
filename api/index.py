"""
Vercel serverless function entry point for FastAPI backend.
This file is used when deploying to Vercel.
"""
from backend.app.main import app

# Vercel expects a handler function
handler = app

