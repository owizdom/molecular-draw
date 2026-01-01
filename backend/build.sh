#!/bin/bash
# Build script for Render.com
# Install RDKit and other requirements

# Install RDKit first
pip install rdkit || pip install --no-cache-dir rdkit-pypi || pip install --no-cache-dir rdkit

# Install other requirements
pip install -r requirements.txt

