# Deployment Guide - Vercel

## Quick Deploy (Easiest Method)

### Option 1: Deploy via Vercel CLI (Recommended)

1. **Install Vercel CLI:**
   ```bash
   npm i -g vercel
   ```

2. **Login to Vercel:**
   ```bash
   vercel login
   ```

3. **Deploy from project root:**
   ```bash
   cd /Users/Apple/Desktop/disambler
   vercel
   ```

4. **Follow the prompts:**
   - Set up and deploy? **Yes**
   - Which scope? (select your account)
   - Link to existing project? **No**
   - Project name? (press Enter for default)
   - Directory? **frontend**
   - Override settings? **No**

5. **Set environment variable:**
   - Go to Vercel Dashboard → Your Project → Settings → Environment Variables
   - Add: `VITE_API_BASE_URL` = `https://your-backend-url.com`

### Option 2: Deploy via GitHub Integration (Automatic)

1. **Go to [vercel.com](https://vercel.com)**
   - Sign up/login with GitHub

2. **Import your repository:**
   - Click "Add New" → "Project"
   - Import `owizdom/molecular-draw`

3. **Configure project:**
   - **Framework Preset:** Vite
   - **Root Directory:** `frontend`
   - **Build Command:** `npm run build`
   - **Output Directory:** `dist`
   - **Install Command:** `npm install`

4. **Add Environment Variable:**
   - In project settings, add:
     - Key: `VITE_API_BASE_URL`
     - Value: Your backend URL (set this after deploying backend)

5. **Deploy:**
   - Click "Deploy"
   - Vercel will automatically deploy on every push to main branch

## Backend Deployment (Free Options)

### Option A: Render.com (Free Tier)

1. Go to [render.com](https://render.com)
2. Sign up with GitHub
3. Click "New" → "Web Service"
4. Connect your repository: `owizdom/molecular-draw`
5. Configure:
   - **Name:** molecular-draw-backend
   - **Root Directory:** `backend`
   - **Environment:** Python 3
   - **Build Command:** `pip install -r requirements.txt`
   - **Start Command:** `uvicorn app.main:app --host 0.0.0.0 --port $PORT`
6. Add environment variables (if needed):
   - `GEMINI_API_KEY` (if using)
   - `HF_API_KEY` (if using)
7. Deploy!

### Option B: Railway.app (Free Tier)

1. Go to [railway.app](https://railway.app)
2. Sign up with GitHub
3. Click "New Project" → "Deploy from GitHub repo"
4. Select your repository
5. Add service → Select `backend` folder
6. Railway auto-detects Python
7. Add environment variables in settings
8. Deploy!

### Option C: Fly.io (Free Tier)

1. Install Fly CLI: `curl -L https://fly.io/install.sh | sh`
2. Login: `fly auth login`
3. In `backend/` directory:
   ```bash
   fly launch
   ```
4. Follow prompts
5. Deploy: `fly deploy`

## After Backend is Deployed

1. **Get your backend URL** (e.g., `https://molecular-draw-backend.onrender.com`)

2. **Update Vercel environment variable:**
   - Go to Vercel Dashboard → Your Project → Settings → Environment Variables
   - Update `VITE_API_BASE_URL` with your backend URL

3. **Redeploy frontend:**
   - Vercel will auto-redeploy, or
   - Go to Deployments → Click "..." → "Redeploy"

## Your App Will Be Live At:
- Frontend: `https://your-project.vercel.app`
- Backend: `https://your-backend-url.com`

## Troubleshooting

- **CORS errors?** Update `backend/app/main.py` to allow your Vercel domain
- **API not working?** Check `VITE_API_BASE_URL` is set correctly
- **Build fails?** Check Vercel build logs for errors

