#!/usr/bin/env bash

# Placeholder SLURM wrapper for Stage 2 of the BRB pipeline.
# This script should call src/brb_pipeline/stage2.py once the implementation is ready.

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <config.yaml>"
    exit 1
fi

CONFIG="$1"

python -m brb_pipeline.stage2 "$CONFIG"
