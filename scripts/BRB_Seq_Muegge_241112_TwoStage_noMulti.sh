#!/bin/bash
#
#SBATCH --job-name=BRB_SEQ
#SBATCH --output=BRB_SEQ_%A_%a.out  # %a inserts the job number
#SBATCH --error=BRB_SEQ_%A_%a.err
#SBATCH --array=1-12%12  # This will create 12 tasks numbered 1-12 and allow 12 concurrent jobs to run. Change 12 to whatever number of jobs you have
#SBATCH --cpus-per-task=4
#SBATCH --mem=75000
#SBATCH --mail-type=END,FAIL

################
# Usage
################

if [ "$#" -lt 3 ] # are there less than three passed command line arguments?                     
then
    echo ""
    echo "error: Too few arguments passed to BRB_Seq_*.sh. You provided $#, 3 required."
    echo ""
    echo "@Usage: sbatch BRB_Seq_Muegge_*.sh PROJECT_DIRECTORY MAPPING_FILE PARAMETERS_FILE"
    echo ""
    echo "     The PROJECT_DIRECTORY will contain all the results of this analysis."
    echo ""
    echo "     The MAPPING_FILE should have a header line,with 5 Tab Separated Columns:"
    echo "           SampleName, Group, RT Barcode Sequence, Index 1 File(s), Index 2 File(s)"
    echo "           - Samplename: name for the sample, only letters, numbers, underscore, or periods"
    echo "           - Group: the main experimental group, e.g. WT or CONTROL. Used for heatmap"
    echo "           - RT Barcode Sequence: the nucleotide sequence of the sample barcode."
    echo "           - Index 1 File(s): path to the file(s) with barcodes. If more than one file,"
    echo "                separate each path by a comma."
    echo "           - Index 2 File(s): path to the file(s) with cDNA sequence for mapping."
    echo "               If more than one file, separate each path by a comma."
    echo ""
    echo "     The PARAMETERS_FILE should be a one-line csv:"
    echo "         Path to Star index, Path to GTF, Path to Genome BED, Demultiplex (Boolean),"
    echo "         Remove intermediate (Boolean)"
    echo "     Note that booleans for Species mixing and Saturation curve analysis have been removed."
    echo ""
    echo "     We recommend that you copy this code to subdirectory scripts in your project directory."
    echo ""
    echo "@Example: sbatch --array=1-8%8 scripts/BRB_Seq_*.sh spike_in_analysis/ metadata/spikein_map.txt metadata/parameters.csv"
    echo ""
    echo "@Example: sbatch --array=1-8%8 scripts/BRB_Seq_*.sh full_run_analysis/ metadata/fullrun_map.txt metadata/parameters.csv"
    echo ""
    exit 1
fi

################
# CHANGE LOG
################

# Adding functions to extract key values from STAR and cutadapt to the log file, to make it easier to summarize the impact of those changes on metrics

# 11/12/2024 - updated fastqc to version 12 in rmlab to address a perl dependency issue caused by htcf system upgrade. Lines 219 and 231

################
# Grab the command line arguments
################

PROJECT_DIR=$1
MAPPING=$2
PARAMETERS=$3

###################
# Capture call parameters and timestamp
###################
now=$(date)
echo 'Calling script ' $0 ' on '$now
echo ""

###################
# Create the project directory if it doesn't already exist
###################
mkdir -p ${PROJECT_DIR}
echo 'Results will be sent to project directory: ' $PROJECT_DIR
echo ""

#############
#  Read-in the data from the mapping file
#############

# The mapping file now includes a header. This is useful in some R scripts. Iterate past it
MAPPING_LINE_READ=$(($SLURM_ARRAY_TASK_ID + 1))

# Extract the 
read samplename group barcode index1 index2 < <( sed -n ${MAPPING_LINE_READ}p $MAPPING )

echo "Processing sample ${samplename}"
echo ""

#############
# Create output directories for the whole project
#############

# The -p option will create the directory only if it does not exist, and it won't throw an error if it does exist.

mkdir -p ${PROJECT_DIR}/
mkdir -p ${PROJECT_DIR}/logfile
mkdir -p ${PROJECT_DIR}/fastqc/
mkdir -p ${PROJECT_DIR}/cutadapt/
mkdir -p ${PROJECT_DIR}/STAR/
mkdir -p ${PROJECT_DIR}/FeatureCounts/
mkdir -p ${PROJECT_DIR}/RSeQC/
mkdir -p ${PROJECT_DIR}/Demultiplexed_Fastq/
mkdir -p ${PROJECT_DIR}/Qualimap/

#############
# Record passed parameters in Log File
#############

LOGFILE=${PROJECT_DIR}/logfile/${samplename}_Log.txt
echo -e "Parameter:\tValue" >$LOGFILE

echo -e "Samplename:\t${samplename}" >>$LOGFILE
echo -e "Group:\t${group}" >>$LOGFILE
echo -e "Sample_barcode:\t${barcode}" >>$LOGFILE
echo -e "Index_1_Files:\t${index1}" >>$LOGFILE
echo -e "Index_2_Files:\t${index2}" >>$LOGFILE

#############
#  Read-in the parameter file
#############

# Need to initialize the variables
STAR_DIR=""
GENOME_GTF=""
GENOME_BED=""
DEMULTIPLEX=""
REMOVE_INTERMEDIATE=""

while IFS="," read -r star_dir genome_gtf genome_bed demultiplex remove_intermediate
do
    { STAR_DIR=$star_dir; }
    { GENOME_GTF=$genome_gtf; }
    { GENOME_BED=$genome_bed; }
    { DEMULTIPLEX=$demultiplex; }
    { REMOVE_INTERMEDIATE=$remove_intermediate; }
done < $PARAMETERS

echo -e "STAR_Directory:\t${STAR_DIR}" >>$LOGFILE
echo -e "Genome_GTF is:\t${GENOME_GTF}" >>$LOGFILE
echo -e "Genome_BED is:\t${GENOME_BED}" >>$LOGFILE
echo -e "Demultiplex_Flag:\t${DEMULTIPLEX}" >>$LOGFILE
echo -e "Remove_Intermediate_Flag:\t${remove_intermediate}" >>$LOGFILE

#############
# Demultiplex by RT Barcode, if Flag is passed
#############

# initialize the read 2 fastq
read2=""

if [ "$DEMULTIPLEX" == "TRUE" ] 
then
    echo "Demultiplexing the index files with Cutadapt:"
    echo ""

    # Note: 2>&1 will redirect stderr (2) to the stdout stream (&1)
    eval $(spack load --sh py-cutadapt@2.10)
    cutadaptVersion="$(cutadapt --version 2>&1)"

    eval $(spack load --sh pigz@2.6)
    pigzVersion="$(pigz --version 2>&1)"
    
    echo -e "Cutadapt_version:\t${cutadaptVersion}" >> $LOGFILE
    echo -e "Cutadapt_Barcode_Demultiplex_Parameters:\t-g=^${barcode} --no-indels --minimum-length 1 --discard-untrimmed -e 0.15 -j ${SLURM_CPUS_PER_TASK}" >>$LOGFILE

    # If there is more than one file in the index, split into an array
    IFS=',' read -r -a array1 <<< "$index1"
    IFS="," read -r -a array2 <<< "$index2"

    for index in "${!array1[@]}"
    do

	filename=$(basename -- ${array2[index]})
	filename="${filename%%.*}"
       	
        echo "${array2[index]}"
	echo "${filename}"
	
	# the adapter is 5' (-g), not 3' (-a)
	# -j will auto-detect number of cores
	cutadapt \
	    -g="^$barcode" \
	    -e 0.15 \
	    --no-indels \
	    --minimum-length 1 \
	    -o /dev/null \
	    -p ${PROJECT_DIR}/cutadapt/${samplename}_${filename}.fq.gz \
	    --discard-untrimmed \
	    --cores ${SLURM_CPUS_PER_TASK} \
	    ${array1[$index]} ${array2[$index]} \
	    >${PROJECT_DIR}/cutadapt/${samplename}_${filename}.log
	
    done
    
    cat ${PROJECT_DIR}/cutadapt/${samplename}_*.fq.gz >${PROJECT_DIR}/Demultiplexed_Fastq/${samplename}_read2.fq.gz
    rm ${PROJECT_DIR}/cutadapt/${samplename}_*.fq.gz
    
    read2=${PROJECT_DIR}/Demultiplexed_Fastq/${samplename}_read2.fq.gz

    eval $(spack unload --sh py-cutadapt@2.10)
    eval $(spack unload --sh pigz@2.6)

else
    read2=$index2
fi

# Generate an md5 value for the demultiplexed file
md5sum ${read2} >${PROJECT_DIR}/Demultiplexed_Fastq/${samplename}_read2.fq.gz.md5sum

#############
# Use FastQC to determine QC metrics
#############

# Note that we only need to operate on Read 2, which has the sequence information
# https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/02_assessing_quality.html

eval $(/ref/rmlab/software/spack-0.22/bin/spack load --sh fastqc@0.12.1)
fastqcVersion="$(fastqc -v 2>&1)"

echo "Analyzing Read 2 input fastq quality with Fastqc."
echo ""

echo -e "Fastqc_version:\t${fastqcVersion}" >>$LOGFILE
echo -e "Fastqc_parameters:\tDefault" >>$LOGFILE

#fastqc only accepts output directory. Samples will be named based on read file
fastqc -t ${SLURM_CPUS_PER_TASK} -o ${PROJECT_DIR}/fastqc/ $read2 2> ${PROJECT_DIR}/fastqc/${samplename}_fastq.sterr.txt

eval $(/ref/rmlab/software/spack-0.22/bin/spack unload --sh fastqc@0.12.1)

#############
# Use CutAdapt to remove short reads and Adapater contaminated sequence
#############

# The adapter sequence passed here is the Illumina Universal Adapter
# Fastq will trim this adapter and anything 3' downstream, which often
#  includes other portions of the library primers

echo "Trimming adapters and short reads from the index file with cutadapt"
echo ""

# Note: 2>&1 will redirect stderr (2) to the stdout stream (&1)
eval $(spack load --sh py-cutadapt@2.10)
cutadaptVersion="$(cutadapt --version 2>&1)"

eval $(spack load --sh pigz@2.6)
pigzVersion="$(pigz --version 2>&1)"

echo -e "Cutadapt_version:\t${cutadaptVersion}" >>$LOGFILE
echo -e "Cutadapt_Adapter_Trim_Parameters:\t--adapter=AGATCGGAAGAG --minimum-length=25 -j ${SLURM_CPUS_PER_TASK}" >>$LOGFILE

cutadapt \
    --adapter=AGATCGGAAGAG \
    --minimum-length=25 \
    -o ${PROJECT_DIR}/cutadapt/${samplename}_temp_trimmed.fq.gz \
    -j ${SLURM_CPUS_PER_TASK} \
    $read2 \
    >${PROJECT_DIR}/cutadapt/${samplename}_AdapterTrim.log

numberDemultReads="$(grep "Total reads processed:" ${PROJECT_DIR}/cutadapt/${samplename}_AdapterTrim.log)"
# extract digits only
numberDemultReads="$(echo ${numberDemultReads} | tr -dc '[:digit:]')"

numberReadsAdapter="$(grep "Reads with adapters:" ${PROJECT_DIR}/cutadapt/${samplename}_AdapterTrim.log)"
# this ends with (##%). Remove that
numberReadsAdapter="$(echo ${numberReadsAdapter} | sed 's/([0-9]\+\.[0-9]\+%)$//g')"
numberReadsAdapter="$(echo ${numberReadsAdapter} | tr -dc '[:digit:]')"

numberReadsWritten="$(grep "Reads written (passing filters):" ${PROJECT_DIR}/cutadapt/${samplename}_AdapterTrim.log)"
# extract digits only
numberReadsWritten="$(echo ${numberReadsWritten} | sed 's/([0-9]\+\.[0-9]\+%)$//g')"
numberReadsWritten="$(echo ${numberReadsWritten} | tr -dc '[:digit:]')"

echo -e "Number demultiplexed reads:\t${numberDemultReads}" >>$LOGFILE
echo -e "Number reads with adapter:\t${numberReadsAdapter}" >>$LOGFILE
echo -e "Number reads written:\t${numberReadsWritten}" >>$LOGFILE

#############
# Use CutAdapt to quantify homopolymer A contamination
#############

# Change Log 1/23/24: we actually remove polyA from sequence here (previously, we were only counting)
# The objective is to see if this improves mapping in Star

cutadapt \
    --adapter="A{30}" \
    --overlap=15 \
    --minimum-length=25 \
    -o ${PROJECT_DIR}/cutadapt/${samplename}_trimmed.fq.gz \
    -j ${SLURM_CPUS_PER_TASK} \
    ${PROJECT_DIR}/cutadapt/${samplename}_temp_trimmed.fq.gz \
    >${PROJECT_DIR}/cutadapt/${samplename}_polyA.log

eval $(spack unload --sh py-cutadapt@2.10)
eval $(spack unload --sh pigz@2.6)

numberReadsPolyA="$(grep "Reads with adapters:" ${PROJECT_DIR}/cutadapt/${samplename}_polyA.log))"
# extract digits only. remove the trailing percentage value
# Note - had to drop the $ anchor here, may not be end of line in polyA report
numberReadsPolyA="$(echo ${numberReadsPolyA} | sed 's/([0-9]\+\.[0-9]\+%)//g')"
numberReadsPolyA="$(echo ${numberReadsPolyA} | tr -dc '[:digit:]')"

echo -e "Number reads with polyA:\t${numberReadsPolyA}" >>$LOGFILE

#############                                                                                                                               
# Use FastQC to determine QC metrics after adapter trimming                                                                                  
#############                                                                                                                                 
# Note that we only need to operate on Read 2, which has the sequence information                                                             
# https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/02_assessing_quality.html                                                  
eval $(/ref/rmlab/software/spack-0.22/bin/spack load --sh fastqc@0.12.1)
fastqcVersion="$(fastqc -v 2>&1)"

echo "Analyzing Read 2 adapter and polyA trimmed sequence quality with Fastqc."
echo ""

echo -e "Fastqc_version:\t${fastqcVersion}" >>$LOGFILE
echo -e "Fastqc_parameters:\tDefault" >>$LOGFILE

#fastqc only accepts output directory. Samples will be named based on read file                                                               
fastqc -t ${SLURM_CPUS_PER_TASK} -o ${PROJECT_DIR}/fastqc/ ${PROJECT_DIR}/cutadapt/${samplename}_trimmed.fq.gz 2> ${PROJECT_DIR}/fastqc/${samplename}_adapterTrimmed_fastq.sterr.txt

eval $(/ref/rmlab/software/spack-0.22/bin/spack unload --sh fastqc@0.12.1)

#############
# Use Star to align the reads
#############

eval $(spack load --sh star@2.7.6a)
starVersion="$(STAR --version 2>&1)"

echo "Aligning reads with STAR"
echo ""

echo -e "STAR_Version:\t${starVersion}" >>$LOGFILE
echo -e "STAR_Parameters:\t--genomeDIR=${STAR_DIR} --outFilterMultiMapNmax 1" >>$LOGFILE

STAR \
    --runThreadN ${SLURM_CPUS_PER_TASK} \
    --genomeDir $STAR_DIR \
    --readFilesCommand zcat \
    --readFilesIn ${PROJECT_DIR}/cutadapt/${samplename}_trimmed.fq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterMultimapNmax 1 \
    --outFileNamePrefix ${PROJECT_DIR}/STAR/${samplename}_

eval $(spack unload --sh star@2.7.6a)

# Pull metrics

numberReadsSTAR="$(grep "Number of input reads" ${PROJECT_DIR}/STAR/${samplename}_Log.final.out)"
# extract digits only
numberReadsSTAR="$(echo ${numberReadsSTAR} | tr -dc '[:digit:]')"

numberUniqueMap="$(grep "Uniquely mapped reads number" ${PROJECT_DIR}/STAR/${samplename}_Log.final.out)"
# extract digits only
numberUniqueMap="$(echo ${numberUniqueMap} | tr -dc '[:digit:]')"

numberTooManyLoci="$(grep "Number of reads mapped to too many loci" ${PROJECT_DIR}/STAR/${samplename}_Log.final.out)"
# extract digits only                                                                                                 
numberTooManyLoci="$(echo ${numberTooManyLoci} | tr -dc '[:digit:]')"

numberTooShort="$(grep "Number of reads unmapped: too short" ${PROJECT_DIR}/STAR/${samplename}_Log.final.out)"
# extract digits only                                                                                                 
numberTooShort="$(echo ${numberTooShort} | tr -dc '[:digit:]')"

echo -e "Number reads input to STAR:\t${numberReadsSTAR}" >>$LOGFILE
echo -e "Number of reads uniquely mapped:\t${numberUniqueMap}" >>$LOGFILE
echo -e "Number of reads too many loci:\t${numberTooManyLoci}" >>$LOGFILE
echo -e "Number of reads too short:\t${numberTooShort}" >>$LOGFILE

#############
# count reads with featureCount from subreads packages
#############

eval $(spack load --sh subread@2.0.2)
featureCountsVersion="$(featureCounts -v 2>&1)"

echo "Running featureCounts using subread"
echo ""

echo -e "FeatureCounts_version:\t${featureCountsVersion}" >>$LOGFILE
echo -e "FeatureCounts_Parameters:\t-a ${GENOME_GTF}" >>$LOGFILE

featureCounts \
   -a $GENOME_GTF \
   -o ${PROJECT_DIR}/FeatureCounts/${samplename}_featureCounts.txt \
   -T ${SLURM_CPUS_PER_TASK} \
   -R BAM \
   ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam \
   2> ${PROJECT_DIR}/FeatureCounts/${samplename}_featurecounts.screen-output.log


# Keep column 1 (gene ID) and the counts data (column 7)
cut -f 1,7 ${PROJECT_DIR}/FeatureCounts/${samplename}_featureCounts.txt \
    >${PROJECT_DIR}/FeatureCounts/${samplename}_featureCounts.RMatrix.txt

# what follows is awesome: sed substitution in place
#https://www.systutorials.com/how-to-delete-a-specific-line-from-a-text-file-in-command-line-on-linux/
#https://thoughtbot.com/blog/sed-102-replace-in-place

# remove the first line which is a header
sed -i '1d' ${PROJECT_DIR}/FeatureCounts/${samplename}_featureCounts.RMatrix.txt

# Replace the first line (contains full path) with samplename
sed -i "1c\Geneid\t${samplename}" ./${PROJECT_DIR}/FeatureCounts/${samplename}_featureCounts.RMatrix.txt

eval $(spack unload --sh subread@2.0.2)

#####################
# post-alignment QC with fastqc
#####################

#echo "Performing post alignment QC with fastqc"

#eval $(spack load --sh fastqc@0.11.9)


#fastqc -t ${SLURM_CPUS_PER_TASK} -o ${PROJECT_DIR}/fastqc/ \
#       ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam \
#       2> ${PROJECT_DIR}/fastqc/${samplename}_Aligned.sortedByCoord.out.bam.sterr.txt

#eval $(spack unload --sh fastqc@0.11.9)

#####################
# QC with RSeqc
#####################

eval $(spack load --sh py-rseqc@2.6.4)
rseqcVersion="$(bam_stat.py -v 2>&1)"
# Hard code for now
#rseqcVersion=${rseqcVersion/bam_stat.py/}
rseqcVersion="2.6.4"

eval $(spack load --sh samtools@1.7)
samtoolsVersion="$(samtools --version 2>&1 | head -1)"


echo "Running QC Analysis with RSeqc"
echo ""

echo -e "Rseqc_version:\t${rsqecVersion}" >>$LOGFILE
echo -e "Samtools_version:\t${samtoolsVersion}" >>$LOGFILE

# Create bai 
samtools index ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam  \
    -@ ${SLURM_CPUS_PER_TASK} \
    ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam.bai \
    2> /dev/null

# Extract alignment information with bam_stats.py
# prints simple file with bam stats. 2 columns with some formatting
bam_stat.py -i ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam \
	    >${PROJECT_DIR}/RSeQC/${samplename}_bam_stats_out.txt \
	    2> /dev/null

# junction saturation. this is relatively fast. Send stdout and stderr to null
junction_saturation.py -i ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam \
		       -r ${GENOME_BED} \
		       -o ${PROJECT_DIR}/RSeQC/${samplename}_junction_saturation \
		       > /dev/null 2> /dev/null

# read_distribution. number of reads in introns, exons, 3' UTR, etc. 
read_distribution.py -i ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam \
		     -r ${GENOME_BED} \
		     > ${PROJECT_DIR}/RSeQC/${samplename}_read_distribution.txt \
		     2> /dev/null
  
# read_duplication
read_duplication.py -i ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam \
		    -o ${PROJECT_DIR}/RSeQC/${samplename}_read_duplication \
		    > /dev/null 2>&1

eval $(spack unload --sh py-rseqc@2.6.4)
eval $(spack unload --sh samtools@1.7)
	   
#####################
# QC with Qualimap
#####################

echo "Running QC Analysis with Qualimap"
echo ""

eval $(spack load --sh qualimap@2.2.1)
# qualimap doesn't have a simple version check, hardcode
qualimapVersion="QualiMap v.2.2.1"

echo -e "Qualimap_version:\t${qualimapVersion}" >>$LOGFILE

#Probably need to dynamically update the ja-mem-size to how much is given by slurm.

qualimap rnaseq \
	 --java-mem-size=70G \
	 -outdir ${PROJECT_DIR}/Qualimap/${samplename} \
	 -gtf $GENOME_GTF \
	 -bam ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam \
	 2> /dev/null

eval $(spack unload --sh qualimap@2.2.1)


#############
# Remove intermediate files
#############

if [ "$REMOVE_INTERMEDIATE" == "TRUE" ] 
then
    # fastqc
	# This file is needed for MultiQC.
    #rm results/fastqc/${samplename}/${samplename}_*fastqc.zip

    # cutadapt
	
    # RSeQC
    #rm 
	
    # STAR
    rm ${PROJECT_DIR}/STAR/${samplename}_SJ.out.tab \
       ${PROJECT_DIR}/STAR/${samplename}_Aligned.sortedByCoord.out.bam \
       ${PROJECT_DIR}/STAR/${samplename}_Log.progress.out \

    # feature count
    rm ${PROJECT_DIR}/FeatureCounts/${samplename}_Aligned.sortedByCoord.out.bam.featureCounts.bam \
       ${PROJECT_DIR}/FeatureCounts/${samplename}_featureCounts.txt
fi
