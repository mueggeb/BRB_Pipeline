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


def trim_adapters(read2_path, sample_name, project_dir, cpus_per_task):
    """
    Trim adapters and short reads from demultiplexed Read 2 using cutadapt.

    Matches the original bash command:
        cutadapt \
            --adapter=AGATCGGAAGAG \
            --minimum-length=25 \
            -o ${PROJECT_DIR}/cutadapt/${samplename}_temp_trimmed.fq.gz \
            -j ${SLURM_CPUS_PER_TASK} \
            $read2 \
            >${PROJECT_DIR}/cutadapt/${samplename}_AdapterTrim.log

    Parameters
    ----------
    read2_path : str or Path
        Path to the demultiplexed Read 2 FASTQ file.
    sample_name : str
        Sample name used for output and log file naming.
    project_dir : str or Path
        Project directory where cutadapt outputs are written.
    cpus_per_task : int
        Number of CPU cores to pass to cutadapt.

    Returns
    -------
    Path
        Path to the adapter-trimmed FASTQ file.

    Raises
    ------
    RuntimeError
        If cutadapt command fails.
    """
    project_dir = Path(project_dir)
    cutadapt_dir = project_dir / "cutadapt"
    cutadapt_dir.mkdir(parents=True, exist_ok=True)

    read2_path = Path(read2_path)
    output_trimmed = cutadapt_dir / f"{sample_name}_temp_trimmed.fq.gz"
    log_file = cutadapt_dir / f"{sample_name}_AdapterTrim.log"

    cmd = [
        "cutadapt",
        "--adapter=AGATCGGAAGAG",
        "--minimum-length=25",
        "-o", str(output_trimmed),
        "-j", str(cpus_per_task),
        str(read2_path),
    ]

    logger.info(f"Running adapter trimming on {read2_path}")
    try:
        with open(log_file, "w") as log_fh:
            subprocess.run(
                cmd,
                stdout=log_fh,
                stderr=subprocess.STDOUT,
                check=True,
                text=True,
            )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Adapter trimming failed for {read2_path}: {e}")

    logger.info(f"Adapter trimming complete, output: {output_trimmed}")
    return output_trimmed


def trim_polyA(read2_path, sample_name, project_dir, cpus_per_task):
    """
    Trim polyA homopolymer sequence from adapter-trimmed Read 2 using cutadapt.

    Matches the original bash command:
        cutadapt \
            --adapter="A{30}" \
            --overlap=15 \
            --minimum-length=25 \
            -o ${PROJECT_DIR}/cutadapt/${samplename}_trimmed.fq.gz \
            -j ${SLURM_CPUS_PER_TASK} \
            ${PROJECT_DIR}/cutadapt/${samplename}_temp_trimmed.fq.gz \
            >${PROJECT_DIR}/cutadapt/${samplename}_polyA.log

    Parameters
    ----------
    read2_path : str or Path
        Path to the adapter-trimmed Read 2 FASTQ file.
    sample_name : str
        Sample name used for output and log file naming.
    project_dir : str or Path
        Project directory where cutadapt outputs are written.
    cpus_per_task : int
        Number of CPU cores to pass to cutadapt.

    Returns
    -------
    Path
        Path to the polyA-trimmed FASTQ file.

    Raises
    ------
    RuntimeError
        If cutadapt command fails.
    """
    project_dir = Path(project_dir)
    cutadapt_dir = project_dir / "cutadapt"
    cutadapt_dir.mkdir(parents=True, exist_ok=True)

    read2_path = Path(read2_path)
    output_trimmed = cutadapt_dir / f"{sample_name}_trimmed.fq.gz"
    log_file = cutadapt_dir / f"{sample_name}_polyA.log"

    cmd = [
        "cutadapt",
        "--adapter=A{30}",
        "--overlap=15",
        "--minimum-length=25",
        "-o", str(output_trimmed),
        "-j", str(cpus_per_task),
        str(read2_path),
    ]

    logger.info(f"Running polyA trimming on {read2_path}")
    try:
        with open(log_file, "w") as log_fh:
            subprocess.run(
                cmd,
                stdout=log_fh,
                stderr=subprocess.STDOUT,
                check=True,
                text=True,
            )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"PolyA trimming failed for {read2_path}: {e}")

    logger.info(f"PolyA trimming complete, output: {output_trimmed}")
    return output_trimmed

def run_star_alignment(trimmed_read2_path, sample_name, project_dir, star_dir, cpus_per_task):
    """
    Align reads with STAR using the adapter- and polyA-trimmed Read 2 FASTQ.

    Matches the original bash command:
        STAR \
            --runThreadN ${SLURM_CPUS_PER_TASK} \
            --genomeDir $STAR_DIR \
            --readFilesCommand zcat \
            --readFilesIn ${PROJECT_DIR}/cutadapt/${samplename}_trimmed.fq.gz \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterMultimapNmax 1 \
            --outFileNamePrefix ${PROJECT_DIR}/STAR/${samplename}_

    Parameters
    ----------
    trimmed_read2_path : str or Path
        Path to the polyA-trimmed Read 2 FASTQ file.
    sample_name : str
        Sample name used for naming STAR output prefix.
    project_dir : str or Path
        Project directory where STAR output is written.
    star_dir : str or Path
        STAR genome directory.
    cpus_per_task : int
        Number of CPU cores to pass to STAR.

    Returns
    -------
    Path
        The STAR output prefix directory path.

    Raises
    ------
    RuntimeError
        If the STAR command fails.
    """
    project_dir = Path(project_dir)
    star_out_dir = project_dir / "STAR"
    star_out_dir.mkdir(parents=True, exist_ok=True)

    trimmed_read2_path = Path(trimmed_read2_path)
    output_prefix = star_out_dir / f"{sample_name}_"

    cmd = [
        "STAR",
        "--runThreadN", str(cpus_per_task),
        "--genomeDir", str(star_dir),
        "--readFilesCommand", "zcat",
        "--readFilesIn", str(trimmed_read2_path),
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outFilterMultimapNmax", "1",
        "--outFileNamePrefix", str(output_prefix),
    ]

    logger.info(f"Running STAR alignment on {trimmed_read2_path}")
    try:
        subprocess.run(
            cmd,
            check=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"STAR alignment failed for {trimmed_read2_path}: {e}")

    logger.info(f"STAR alignment complete, output prefix: {output_prefix}")
    return output_prefix


def run_featurecounts(star_output_prefix, sample_name, project_dir, genome_gtf, cpus_per_task):
    """
    Run featureCounts on the STAR-aligned BAM file.

    Matches the original bash command:
        featureCounts \
           -a $GENOME_GTF \
           -o ${PROJECT_DIR}/FeatureCounts/${samplename}_featureCounts.txt \
           -T ${SLURM_CPUS_PER_TASK} \
           -R BAM \
           ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam \
           2> ${PROJECT_DIR}/FeatureCounts/${samplename}_featurecounts.screen-output.log

    Parameters
    ----------
    star_output_prefix : str or Path
        STAR output prefix path returned by run_star_alignment.
    sample_name : str
        Sample name used for output and log file naming.
    project_dir : str or Path
        Project directory where featureCounts output is written.
    genome_gtf : str or Path
        Path to the GTF annotation file.
    cpus_per_task : int
        Number of CPU cores to pass to featureCounts.

    Returns
    -------
    Path
        Path to the featureCounts RMatrix output file.

    Raises
    ------
    RuntimeError
        If the featureCounts command fails.
    """
    project_dir = Path(project_dir)
    featurecounts_dir = project_dir / "FeatureCounts"
    featurecounts_dir.mkdir(parents=True, exist_ok=True)

    star_output_prefix = Path(star_output_prefix)
    bam_path = project_dir / "STAR" / f"{sample_name}_Aligned.sortedByCoord.out.bam"
    featurecounts_txt = featurecounts_dir / f"{sample_name}_featureCounts.txt"
    screen_output_log = featurecounts_dir / f"{sample_name}_featurecounts.screen-output.log"
    rmatrix_path = featurecounts_dir / f"{sample_name}_featureCounts.RMatrix.txt"

    cmd = [
        "featureCounts",
        "-a", str(genome_gtf),
        "-o", str(featurecounts_txt),
        "-T", str(cpus_per_task),
        "-R", "BAM",
        str(bam_path),
    ]

    logger.info(f"Running featureCounts on {bam_path}")
    try:
        with open(screen_output_log, "w") as stderr_fh:
            subprocess.run(
                cmd,
                stderr=stderr_fh,
                check=True,
                text=True,
            )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"featureCounts failed for {bam_path}: {e}")

    logger.info(f"featureCounts complete, output: {featurecounts_txt}")
    return rmatrix_path


def run_fastqc_post_alignment(sample_name, project_dir, cpus_per_task):
    """
    Run FastQC on the STAR-aligned BAM file for post-alignment QC.

    Matches the original bash command:
        fastqc -t ${SLURM_CPUS_PER_TASK} -o ${PROJECT_DIR}/fastqc/ \
               ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam \
               2> ${PROJECT_DIR}/fastqc/${samplename}_Aligned.sortedByCoord.out.bam.sterr.txt

    Parameters
    ----------
    sample_name : str
        Sample name used for output and stderr log naming.
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
        If FastQC command fails.
    """
    project_dir = Path(project_dir)
    fastqc_dir = project_dir / "fastqc"
    fastqc_dir.mkdir(parents=True, exist_ok=True)

    bam_path = project_dir / "STAR" / f"{sample_name}_Aligned.sortedByCoord.out.bam"
    stderr_path = fastqc_dir / f"{sample_name}_Aligned.sortedByCoord.out.bam.sterr.txt"

    cmd = [
        "fastqc",
        "-t", str(cpus_per_task),
        "-o", str(fastqc_dir),
        str(bam_path),
    ]

    logger.info(f"Running post-alignment FastQC on {bam_path}")
    try:
        with open(stderr_path, "w") as stderr_fh:
            subprocess.run(
                cmd,
                stderr=stderr_fh,
                check=True,
                text=True,
            )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Post-alignment FastQC failed for {bam_path}: {e}")

    logger.info(f"Post-alignment FastQC complete, output directory: {fastqc_dir}")
    return fastqc_dir


def main(config_path, sample_index):
    """Run stage 1 for a single sample."""
    raise NotImplementedError("stage1 logic not implemented yet")
