"""Stage 1: per-sample BRB-seq processing."""

import logging
import subprocess
import tempfile
from pathlib import Path


logger = logging.getLogger(__name__)


def demultiplex_by_barcode(barcode, sample_name, index1_paths, index2_paths, 
                           project_dir, cpus_per_task):
    """
    Demultiplex BRB-seq reads by RT barcode using cutadapt.

    Runs cutadapt on paired-end reads, discarding Read 1 (barcode) and keeping
    only Read 2 (cDNA) that match the provided barcode sequence at the 5' end.
    Multiple index files are supported (comma-separated) and concatenated.

    Parameters
    ----------
    barcode : str
        The RT barcode sequence to match at 5' end (e.g., "ATCG").
    sample_name : str
        Sample name for naming output files.
    index1_paths : str
        Path(s) to Read 1 file(s), comma-separated if multiple.
    index2_paths : str
        Path(s) to Read 2 file(s), comma-separated if multiple.
    project_dir : Path or str
        Project directory where output will be written.
    cpus_per_task : int
        Number of CPU cores to use for cutadapt.

    Returns
    -------
    Path
        Path to demultiplexed Read 2 fastq.gz file.

    Raises
    ------
    RuntimeError
        If cutadapt command fails.
    """
    project_dir = Path(project_dir)
    cutadapt_dir = project_dir / "cutadapt"
    demultiplexed_dir = project_dir / "Demultiplexed_Fastq"
    
    cutadapt_dir.mkdir(parents=True, exist_ok=True)
    demultiplexed_dir.mkdir(parents=True, exist_ok=True)

    # Parse comma-separated file lists into arrays
    index1_list = [f.strip() for f in index1_paths.split(",")]
    index2_list = [f.strip() for f in index2_paths.split(",")]

    if len(index1_list) != len(index2_list):
        raise ValueError(
            f"Number of Read 1 files ({len(index1_list)}) "
            f"does not match Read 2 files ({len(index2_list)})"
        )

    cutadapt_outputs = []

    # Process each paired file
    for index1_file, index2_file in zip(index1_list, index2_list):
        index1_path = Path(index1_file)
        index2_path = Path(index2_file)

        # Extract base filename without extension for output naming
        filename_base = index2_path.stem
        if filename_base.endswith(".fastq"):
            filename_base = filename_base[:-6]

        output_r2 = cutadapt_dir / f"{sample_name}_{filename_base}.fq.gz"
        log_file = cutadapt_dir / f"{sample_name}_{filename_base}.log"

        logger.info(f"Demultiplexing with barcode {barcode}: {index2_path}")

        # Run cutadapt with exact parameters from bash script
        cmd = [
            "cutadapt",
            f"-g=^{barcode}",
            "-e", "0.15",
            "--no-indels",
            "--minimum-length", "1",
            "-o", "/dev/null",
            "-p", str(output_r2),
            "--discard-untrimmed",
            "--cores", str(cpus_per_task),
            str(index1_path),
            str(index2_path),
        ]

        try:
            with open(log_file, "w") as log_fh:
                result = subprocess.run(
                    cmd,
                    stdout=log_fh,
                    stderr=subprocess.STDOUT,
                    check=True,
                    text=True,
                )
            logger.info(f"Cutadapt completed for {index2_path}, log: {log_file}")
            cutadapt_outputs.append(output_r2)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"cutadapt failed for {index2_path}: {e}"
            )

    # Concatenate all demultiplexed outputs
    final_output = demultiplexed_dir / f"{sample_name}_read2.fq.gz"
    logger.info(f"Concatenating demultiplexed files into {final_output}")

    concat_cmd = ["cat"] + [str(f) for f in cutadapt_outputs]
    with open(final_output, "wb") as out_fh:
        subprocess.run(
            concat_cmd,
            stdout=out_fh,
            check=True,
        )

    # Clean up intermediate demultiplexed files
    for intermediate in cutadapt_outputs:
        intermediate.unlink()
        logger.info(f"Removed intermediate: {intermediate}")

    # Compute MD5 checksum
    md5_file = str(final_output) + ".md5sum"
    md5_cmd = ["md5sum", str(final_output)]
    with open(md5_file, "w") as md5_fh:
        subprocess.run(md5_cmd, stdout=md5_fh, check=True)

    logger.info(f"Demultiplexing complete. Output: {final_output}")
    return final_output


def run_fastqc(read2_path, sample_name, project_dir, cpus_per_task):
    """
    Run FastQC on Read 2 input.

    Matches the original bash command:
        fastqc -t ${SLURM_CPUS_PER_TASK} -o ${PROJECT_DIR}/fastqc/ $read2

    Parameters
    ----------
    read2_path : str or Path
        Path to the demultiplexed Read 2 FASTQ file.
    sample_name : str
        Sample name for naming stderr log file.
    project_dir : str or Path
        Project directory where FastQC output is written.
    cpus_per_task : int
        Number of CPU cores to pass to FastQC.

    Returns
    -------
    Path
        Path to the FastQC output directory.

    Raises
    ------
    RuntimeError
        If the FastQC command fails.
    """
    project_dir = Path(project_dir)
    fastqc_dir = project_dir / "fastqc"
    fastqc_dir.mkdir(parents=True, exist_ok=True)

    read2_path = Path(read2_path)
    stderr_path = fastqc_dir / f"{sample_name}_fastq.sterr.txt"

    cmd = [
        "fastqc",
        "-t", str(cpus_per_task),
        "-o", str(fastqc_dir),
        str(read2_path),
    ]

    logger.info(f"Running FastQC on {read2_path} with {cpus_per_task} threads")

    try:
        with open(stderr_path, "w") as stderr_fh:
            subprocess.run(
                cmd,
                stderr=stderr_fh,
                check=True,
                text=True,
            )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"FastQC failed for {read2_path}: {e}")

    logger.info(f"FastQC complete, output directory: {fastqc_dir}")
    return fastqc_dir


def main(config_path, sample_index):
    """Run stage 1 for a single sample."""
    raise NotImplementedError("stage1 logic not implemented yet")
