# GitHub Pages Deployment Guide

## Quick Setup

1. **Enable GitHub Pages in your repository:**
   - Go to your repository: https://github.com/owizdom/molecular-draw
   - Click **Settings** → **Pages**
   - Under **Source**, select **GitHub Actions**
   - Save

2. **The deployment will happen automatically:**
   - Every push to `main` branch will trigger a deployment
   - Or manually trigger via **Actions** tab → **Deploy to GitHub Pages** → **Run workflow**

3. **Your site will be available at:**
   - `https://owizdom.github.io/molecular-draw/`

## Configuration

The frontend is already configured to:
- Use the Render backend URL (`https://molecular-draw.onrender.com`) in production
- Set the correct base path (`/molecular-draw/`) for GitHub Pages
- Automatically deploy on every push to main

## Manual Deployment

If you want to deploy manually:

```bash
cd frontend
npm install
GITHUB_PAGES=true VITE_API_BASE_URL=https://molecular-draw.onrender.com npm run build
```

## Troubleshooting

- **404 errors**: Make sure GitHub Pages is set to use **GitHub Actions** as the source
- **API not working**: Check that the Render backend is running at `https://molecular-draw.onrender.com`
- **Build fails**: Check the Actions tab for error logs

