#!/usr/bin/env bash

# Placeholder SLURM wrapper for Stage 1 of the BRB pipeline.
# This script should call src/brb_pipeline/stage1.py once the implementation is ready.

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <config.yaml> <sample-index>"
    exit 1
fi

CONFIG="$1"
SAMPLE_INDEX="$2"

python -m brb_pipeline.stage1 "$CONFIG" "$SAMPLE_INDEX"
