#!/usr/bin/env bash
#SBATCH --job-name=BRB_SEQ
#SBATCH --output=BRB_SEQ_%A_%a.out
#SBATCH --error=BRB_SEQ_%A_%a.err
#SBATCH --array=1-12%12
#SBATCH --cpus-per-task=4
#SBATCH --mem=75000
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config.yaml>"
    exit 1
fi

CONFIG="$1"

# Activate the shared mamba environment for the BRB pipeline.
# You must initialize mamba in your .bashrc for this to work.
eval "$(mamba shell hook --shell bash)"
mamba activate /ref/bmlab/software/envs/brb_pipeline

PIPELINE_DIR=/ref/bmlab/software/brb_python_2026/BRB_Pipeline
python "$PIPELINE_DIR/src/brb_pipeline/stage1.py" "$CONFIG" "$SLURM_ARRAY_TASK_ID"
