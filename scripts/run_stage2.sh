#!/usr/bin/env bash
#SBATCH --job-name=BRB_SEQ_STAGE2
#SBATCH --output=BRB_SEQ_Stage2_%A.out
#SBATCH --error=BRB_SEQ_Stage2_%A.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=75000
#SBATCH --mail-type=END,FAIL

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config.yaml>"
    echo "Note: Run this script only after all Stage 1 SLURM array jobs have completed successfully."
    exit 1
fi

CONFIG="$1"

# Activate the shared mamba environment for the BRB pipeline.
eval "$(mamba shell hook --shell bash)"
mamba activate /ref/bmlab/software/envs/brb_pipeline

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python "$PIPELINE_DIR/src/brb_pipeline/stage2.py" "$CONFIG"
