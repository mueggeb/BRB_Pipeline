"""Stage 2: aggregate BRB-seq results and run MultiQC."""

import logging
import subprocess
import shutil
import argparse
import sys
import pandas as pd
from pathlib import Path


logger = logging.getLogger(__name__)

# Module-level flag for dry-run mode
_dry_run = False


def _run_command(cmd, stdout=None, stderr=None, **kwargs):
    """
    Execute a command or print it in dry-run mode.

    stdout and stderr accept Path objects and will be opened automatically,
    since subprocess.run requires open file objects rather than paths.

    Parameters
    ----------
    cmd : list
        Command as list of strings (for subprocess.run).
    stdout : Path, str, or special constant, optional
        Where to send stdout. Path/str values are opened as files.
    stderr : Path, str, or special constant, optional
        Where to send stderr. Path/str values are opened as files.
    **kwargs
        Additional arguments to pass to subprocess.run (e.g. check, text).
    """
    global _dry_run
    if _dry_run:
        logger.info(f"[DRY RUN] {' '.join(str(c) for c in cmd)}")
        return

    file_mode = 'w' if kwargs.get('text', False) else 'wb'
    open_files = []
    try:
        if isinstance(stdout, (str, Path)):
            f = open(stdout, file_mode)
            open_files.append(f)
            stdout = f
        if isinstance(stderr, (str, Path)):
            f = open(stderr, file_mode)
            open_files.append(f)
            stderr = f
        subprocess.run(cmd, stdout=stdout, stderr=stderr, **kwargs)
    finally:
        for f in open_files:
            f.close()


def _makedirs(path):
    """
    Create directory or print in dry-run mode.

    Parameters
    ----------
    path : Path or str
        Directory path to create.
    """
    global _dry_run
    path = Path(path)
    if _dry_run:
        logger.info(f"[DRY RUN] mkdir -p {path}")
        return
    path.mkdir(parents=True, exist_ok=True)


def _write_file(path, content):
    """
    Write content to file or print in dry-run mode.

    Parameters
    ----------
    path : Path or str
        File path to write.
    content : str
        Content to write.
    """
    global _dry_run
    path = Path(path)
    if _dry_run:
        logger.info(f"[DRY RUN] write {len(content)} bytes to {path}")
        return
    with open(path, "w") as fh:
        fh.write(content)


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
    _makedirs(multiqc_dir)

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
            template_text = template_path.read_text()
            rendered_template = template_text.replace("${LIBRARY_NAME}", library_name)
            rendered_config_path = multiqc_dir / "multiqc_config_rendered.yaml"
            _write_file(rendered_config_path, rendered_template)
            cmd.extend(["--config", str(rendered_config_path)])

    log_file = multiqc_dir / "multiqc.log"
    logger.info(f"Running MultiQC with outputs in {multiqc_dir}")

    try:
        _run_command(
            cmd,
            stderr=log_file,
            stdout=log_file,
            check=True,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"MultiQC failed: {exc}")

    report_path = multiqc_dir / f"{library_name}_QCReport.html"
    logger.info(f"MultiQC generated report at {report_path}")
    return report_path


def cleanup_stage2(project_dir, remove_intermediate):
    """
    Remove per-sample intermediate outputs produced in stage 1.

    This runs only when the remove_intermediate flag is True in the config.
    Files removed match the original bash cleanup behavior.
    """
    if not remove_intermediate:
        logger.info("Skipping stage 2 cleanup because remove_intermediate is False")
        return

    project_dir = Path(project_dir)
    cleanup_patterns = [
        "RSeQC/*",
        "fastqc/*",
        "cutadapt/*.fq.gz",
        "Demultiplexed_Fastq/*",
        "STAR/*Log.out",
        "STAR/*bam.bai",
        "FeatureCounts/*RMatrix.txt",
        "FeatureCounts/*screen-output.log",
    ]

    for pattern in cleanup_patterns:
        for path in project_dir.glob(pattern):
            if _dry_run:
                logger.info(f"[DRY RUN] rm {path}")
                continue
            try:
                if path.is_dir():
                    shutil.rmtree(path)
                    logger.info(f"Removed directory: {path}")
                else:
                    path.unlink()
                    logger.info(f"Removed file: {path}")
            except OSError as exc:
                logger.warning(f"Could not remove {path}: {exc}")


def validate_stage1_outputs(project_dir):
    """
    Preflight validation: check for expected Stage 1 output files.

    Ensures that Stage 1 completed successfully by checking for the presence
    of key output files. Raises a clear error if any are missing.

    Parameters
    ----------
    project_dir : Path or str
        Project directory to validate.

    Raises
    ------
    RuntimeError
        If expected Stage 1 outputs are missing.
    """
    project_dir = Path(project_dir)

    # Check for at least one RMatrix file
    rmatrix_files = list(project_dir.glob("FeatureCounts/*_featureCounts.RMatrix.txt"))
    if not rmatrix_files:
        raise RuntimeError(
            f"No RMatrix files found in {project_dir}/FeatureCounts/. "
            "Check that Stage 1 completed successfully."
        )

    # Check for at least one STAR Log.final.out
    star_logs = list(project_dir.glob("STAR/*Log.final.out"))
    if not star_logs:
        raise RuntimeError(
            f"No STAR alignment logs found in {project_dir}/STAR/. "
            "Check that Stage 1 completed successfully."
        )

    # Check for at least one cutadapt log
    cutadapt_logs = list(project_dir.glob("cutadapt/*.log"))
    if not cutadapt_logs:
        raise RuntimeError(
            f"No cutadapt logs found in {project_dir}/cutadapt/. "
            "Check that Stage 1 completed successfully."
        )

    logger.info(
        f"Stage 1 preflight validation passed: {len(rmatrix_files)} samples, "
        f"{len(star_logs)} STAR logs, {len(cutadapt_logs)} cutadapt logs found."
    )


def main(config_path, dry_run=False):
    """Run stage 2 aggregation and reporting."""
    global _dry_run
    _dry_run = dry_run

    from brb_pipeline.config import load_config

    logger.info(f"Loading config from {config_path}")
    config = load_config(config_path)
    project_dir = config.output_dir
    library_name = config.library_name

    logger.info(f"Starting Stage 2 for project {config.project_name}")
    validate_stage1_outputs(project_dir)
    merge_featurecounts(project_dir, library_name)
    run_multiqc(project_dir, library_name, config.multiqc_template)
    cleanup_stage2(project_dir, config.remove_intermediate)

    logger.info("Stage 2 completed successfully")


if __name__ == "__main__":
    # Add src directory to Python path for absolute imports
    src_dir = Path(__file__).parent.parent
    sys.path.insert(0, str(src_dir))

    parser = argparse.ArgumentParser(
        description="Stage 2: aggregate BRB-seq results and run MultiQC."
    )
    parser.add_argument("config", help="Path to YAML config file")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands and file operations without executing them",
    )
    args = parser.parse_args()
    main(args.config, dry_run=args.dry_run)
