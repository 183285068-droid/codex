#!/usr/bin/env bash
set -eu
cd /workspace/codex
mkdir -p /workspace/scratch/matplotlib /workspace/scratch/numba-cache /workspace/scratch/jupyter /workspace/scratch/jupyter-runtime /workspace/scratch/ipython /workspace/scratch/pip-cache
if [ ! -x /workspace/scratch/pgam5-venv/bin/python ]; then
  python -m venv --system-site-packages /workspace/scratch/pgam5-venv
fi
export MPLCONFIGDIR=/workspace/scratch/matplotlib
export NUMBA_CACHE_DIR=/workspace/scratch/numba-cache
export PIP_CACHE_DIR=/workspace/scratch/pip-cache
/workspace/scratch/pgam5-venv/bin/python -m pip install -r results/GSE202642_HCC_allcell_annotation/requirements.txt
export JUPYTER_DATA_DIR=/workspace/scratch/jupyter
export IPYTHONDIR=/workspace/scratch/ipython
/workspace/scratch/pgam5-venv/bin/python -m ipykernel install --user --name pgam5 --display-name 'Python (PGAM5)'
/workspace/scratch/pgam5-venv/bin/python -c 'import scanpy, anndata, igraph, harmonypy; print("Single-cell runtime imports passed")'
