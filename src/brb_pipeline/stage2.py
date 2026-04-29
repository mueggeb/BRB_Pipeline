"""Stage 2: aggregate BRB-seq results and run MultiQC."""

import logging
import pandas as pd
from pathlib import Path


logger = logging.getLogger(__name__)


def merge_featurecounts(project_dir, library_name, sequencing_type="Full_Run"):
    """
    Merge per-sample FeatureCounts RMatrix files into a combined count matrix.

    Uses pandas to read and merge the per-sample RMatrix files on the Geneid column.
    Each RMatrix file contains Geneid and count columns for one sample.
    """
    project_dir = Path(project_dir)
    featurecounts_dir = project_dir / "FeatureCounts"
    output_file = project_dir / f"{library_name}_{sequencing_type}_FeatureCounts.txt"

    rmatrix_files = list(featurecounts_dir.glob("*_featureCounts.RMatrix.txt"))
    if not rmatrix_files:
        raise RuntimeError(f"No RMatrix files found in {featurecounts_dir}")

    logger.info(f"Merging {len(rmatrix_files)} FeatureCounts RMatrix files")

    # Read and merge all RMatrix files on Geneid column
    dfs = []
    for rmatrix_file in rmatrix_files:
        df = pd.read_csv(rmatrix_file, sep='\t')
        dfs.append(df)

    # featureCounts should guarantee identical gene order for files generated
    # using the same GTF. We still use an outer merge defensively in case
    # different GTF versions were used across samples or row order varies.
    if len(dfs) == 1:
        merged_df = dfs[0]
    else:
        merged_df = dfs[0]
        for df in dfs[1:]:
            merged_df = pd.merge(merged_df, df, on='Geneid', how='outer')

    # Fill NaN values with 0 (genes not present in some samples)
    merged_df = merged_df.fillna(0)

    # Convert count columns to integers
    count_columns = [col for col in merged_df.columns if col != 'Geneid']
    merged_df[count_columns] = merged_df[count_columns].astype(int)

    # Write the merged dataframe to output file
    merged_df.to_csv(output_file, sep='\t', index=False)

    logger.info(f"FeatureCounts merged into {output_file} with {len(merged_df)} genes and {len(count_columns)} samples")
    return output_file


def run_multiqc(project_dir, library_name, config_template_path=None):
    """
    Run MultiQC on per-sample QC and log outputs.

    Passes the expected directories/files to MultiQC and exports aggregated
    metrics as TSV flat files via --data-format tsv.
    """
    project_dir = Path(project_dir)
    multiqc_dir = project_dir / "MultiQC"
    multiqc_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "multiqc",
        str(project_dir / "fastqc" / "*fastqc.zip"),
        str(project_dir / "cutadapt" / "*.log"),
        str(project_dir / "STAR" / "*Log.final.out"),
        str(project_dir / "FeatureCounts" / "*featureCounts.txt.summary"),
        str(project_dir / "RSeQC" / "*"),
        "--data-format", "tsv",
        "--filename", f"{library_name}_QCReport",
        "--outdir", str(multiqc_dir),
        "-f",
    ]

    if config_template_path:
        template_path = Path(config_template_path)
        if template_path.exists():
            cmd.extend(["--config", str(template_path)])

    log_file = multiqc_dir / "multiqc.log"
    logger.info(f"Running MultiQC with outputs in {multiqc_dir}")

    try:
        with open(log_file, "w") as log_fh:
            subprocess.run(cmd, stderr=log_fh, stdout=log_fh, check=True, text=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"MultiQC failed: {exc}")

    report_path = multiqc_dir / f"{library_name}_QCReport.html"
    logger.info(f"MultiQC generated report at {report_path}")
    return report_path


def main(config_path):
    """Run stage 2 aggregation and reporting."""
    raise NotImplementedError("stage2 logic not implemented yet")
