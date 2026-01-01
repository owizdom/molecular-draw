# Molecular Draw API - Backend

FastAPI backend with Celery for asynchronous task processing.

## Setup

1. Install dependencies:
```bash
pip install -r requirements.txt
```

2. Set up Redis (message broker for Celery):
```bash
# macOS
brew install redis
brew services start redis

# Linux
sudo apt-get install redis-server
sudo systemctl start redis

# Or use Docker
docker run -d -p 6379:6379 redis:alpine
```

3. Set up environment variables:
```bash
cp .env.example .env
# Edit .env and add:
# OPENAI_API_KEY=your_key_here
# CELERY_BROKER_URL=redis://localhost:6379/0 (optional, defaults to this)
# CELERY_RESULT_BACKEND=redis://localhost:6379/0 (optional, defaults to this)
```

4. Run the FastAPI server:
```bash
uvicorn app.main:app --reload
```

5. Run the Celery worker (in a separate terminal):
```bash
celery -A app.celery_app worker --loglevel=info
```

The API will be available at `http://localhost:8000`

## API Endpoints

### Synchronous Endpoints
- `POST /structure/from-smiles?async_task=false` - Convert SMILES to 3D structure (immediate)
- `POST /structure/to-smiles` - Convert structure to SMILES
- `POST /structure/properties?async_task=false` - Get molecular properties (immediate)
- `POST /structure/validate?async_task=false` - Validate molecular structure (immediate)
- `WebSocket /ws/chat` - Real-time chat via WebSocket

### Asynchronous Endpoints (Celery)
- `POST /structure/from-smiles?async_task=true` - Convert SMILES to 3D structure (async)
- `POST /chat?async_task=true` - Chat with AI assistant (async, default)
- `POST /structure/generate?async_task=true` - Generate structure from description (async, default)
- `POST /structure/optimize?async_task=true` - Optimize molecular structure (async, default)
- `POST /structure/properties?async_task=true` - Get molecular properties (async)
- `POST /structure/validate?async_task=true` - Validate molecular structure (async)
- `GET /tasks/{task_id}` - Get task status and result

### Task Status Response
When using async endpoints, you'll receive:
```json
{
  "task_id": "abc123...",
  "status": "processing",
  "message": "Task submitted"
}
```

Poll the task status:
```bash
GET /tasks/{task_id}
```

Response when complete:
```json
{
  "task_id": "abc123...",
  "status": "completed",
  "state": "SUCCESS",
  "result": { ... }
}
```

## Dependencies

- FastAPI - Web framework
- RDKit - Molecular processing
- OpenAI - AI assistant (optional, set OPENAI_API_KEY)

