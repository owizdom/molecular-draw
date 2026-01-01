# Molecular Draw

AI-powered molecular modeling application with 3D visualization and natural language interface.

## Features

- 3D molecular visualization with React Three Fiber
- AI-powered natural language interface (OpenAI, Gemini, Hugging Face)
- Structure editor with 118 periodic table elements
- Real-time task progress with Server-Sent Events
- Export functionality (MOL, SDF, JSON)
- 200+ common molecules in instant lookup dictionary
- Asynchronous task processing with Celery and Redis

## Prerequisites

- Python 3.9+
- Node.js 18+
- Redis server
- Python virtual environment

## Installation

### Backend

1. Navigate to the backend directory:
```bash
cd backend
```

2. Create and activate a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

4. Set up environment variables:
Create a `.env` file in the backend directory with:
```
OPENAI_API_KEY=your_openai_key_here
GEMINI_API_KEY=your_gemini_key_here
HUGGINGFACE_API_KEY=your_huggingface_key_here
```

### Frontend

1. Navigate to the frontend directory:
```bash
cd frontend
```

2. Install dependencies:
```bash
npm install
```

## Running the Application

### Start Redis

Make sure Redis is running:
```bash
redis-server
```

### Start Backend

1. Activate the virtual environment:
```bash
cd backend
source venv/bin/activate
```

2. Start the Celery worker:
```bash
celery -A app.celery_app worker --loglevel=info
```

3. In a separate terminal, start the FastAPI server:
```bash
cd backend
source venv/bin/activate
python -m uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```

### Start Frontend

1. Navigate to the frontend directory:
```bash
cd frontend
```

2. Start the development server:
```bash
npm run dev
```

The application will be available at `http://localhost:5173` (or the port shown in the terminal).

## Usage

1. Open the application in your browser
2. Use the AI chat to create molecules by typing natural language commands like:
   - "Create methane"
   - "Create aspirin"
   - "Create benzene"
3. View the 3D structure in the molecular viewer
4. Edit structures using the structure editor
5. Export structures in MOL, SDF, or JSON format
