#!/bin/bash
#
#SBATCH --job-name=Multiplex_BRB
#SBATCH --output=Multiplex_BRB_%A_%a.out  # %a inserts the job number
#SBATCH --error=Multiplex_BRB_%A_%a.err
#SBATCH --array=1%1  
#SBATCH --cpus-per-task=8
#SBATCH --mem=75000
#SBATCH --mail-type=END,FAIL

################
# Usage
################

if [ "$#" -lt 7 ] # are there less than three passed command line arguments?                     
then
    echo ""
    echo "Purpose: This script takes multiplexed BRB-seq style fastq files,"
    echo "  aligns reads using StarSolo, and generates a count matrix for"
    echo "  downstream RNA-seq analysis."
    echo ""
    echo "ERROR: Too few arguments passed to Multiplex_BRB_Align_Counts.sh. You provided $#, 7 required."
    echo ""
    echo "@Usage: sbatch  Multiplex_BRB_Align_Counts.sh PROJECT_DIRECTORY PROJECT_STEM MAPPING_FILE PARAMETERS SPECIES READ1_FASTQ READ2_FASTQ"
    echo ""
    echo "     The PROJECT_DIRECTORY will contain all the results of this analysis."
    echo ""
    echo "     PROJECT_STEM is a short string to uniquely label the project, e.g. BDM07"
    echo ""
    echo "     The MAPPING_FILE should have a header line with 2 Tab Separated Columns:"
    echo "           SampleName, RT Barcode Sequence"
    echo "           - Samplename: name for the sample, only letters, numbers, underscore, or periods"
    echo "           - RT Barcode Sequence: the nucleotide sequence of the sample barcode."
    echo ""
    echo "     The Parameters file should have a header line with 4 Tab Separated Columns:"
    echo "           Cell Barcode Start Position, Cell Barcode Length, UMI Start Position, UMI Length"
    echo ""
    echo "           - Cell_BC and UMI value denotes start position and length in read 1."
    echo "               Prime-Seq: CB_Start=1, CB_Len=12, UMI_Start=13, UMI_Len=16"
    echo "               Alithea  : CB_Start=1, CB_Len=14, UMI_Start=15, UMI_Len=14"
    echo "               Tripp BRB: CB_Start=1, CB_Len=16, UMI_Start=17, UMI_Len=10"
    echo ""
    echo "     Species must be exactly human or mouse. These are hardcoded in the star script."
    echo ""
    echo "     READ1_FASTQ and READ2_FASTQ are the files returned from the sequencing core."
    echo "          - Read1 contains Cell Barcode and UMI"
    echo "          - Read2 contains sequence to align containing data"
    echo ""
    echo "     We recommend that you copy this script to a subdirectory called \"scripts\" in your project directory."
    echo ""
    echo "@Example: sbatch scripts/Multiplex_BRB_Align_Counts.sh FullRun_Analysis/ BDM01_Stem metadata/FullRun_Mapping.txt Alithea_Parameters.txt mouse <path to read 1> <path to read 2>"
    echo ""
    exit 1
fi

################
# CHANGE LOG
################

# Send output of R script directly to log file instead of standard out
# add two stage adapter trimming before star solo, as we found polyA trimming significantly increases number of mapped reads

################
# Grab the command line arguments
################

PROJECT_DIR=$1
PROJECT_STEM=$2
MAPPING=$3
PARAMETERS=$4
SPECIES=$5
READ1=$6
READ2=$7

# TURN THIS OFF FOR DEBUGGING
REMOVE_BAM_FLAG="TRUE"

###################
# Capture call parameters and timestamp
###################
now=$(date)
echo 'Calling script ' $0 ' on '$now
echo ""

###################
# Create the project directory if it doesn't already exist
###################

# Use shell parameter expansion to clear any trailing slashes
PROJECT_DIR="${PROJECT_DIR%/}"

mkdir -p ${PROJECT_DIR}
echo 'Results will be sent to project directory: ' $PROJECT_DIR
echo ""

#############
# Create output directories for the whole project
#############

# The -p option will create the directory only if it does not exist, and it won't throw an error if it does exist.

mkdir -p ${PROJECT_DIR}/Star/
mkdir -p ${PROJECT_DIR}/cutadapt/
mkdir -p ${PROJECT_DIR}/Counts_Files/

#############
# Make Logfile for recording run parameters and results
#############

LOGFILE=${PROJECT_DIR}/${PROJECT_STEM}_Log.txt
echo -e "Parameter:\tValue" >$LOGFILE

echo -e "Project_Directory:\t${PROJECT_DIR}" >>$LOGFILE
echo -e "Project_Stem:\t${PROJECT_STEM}" >>$LOGFILE
echo -e "Mapping_File:\t${MAPPING}" >>$LOGFILE
echo -e "Parameters_File:\t${PARAMETERS}" >>$LOGFILE
echo -e "Species:\t${SPECIES}" >>$LOGFILE
echo -e "Read1_File:\t${READ1}" >>$LOGFILE
echo -e "Read2_File:\t${READ2}" >>$LOGFILE
echo -e "Remove_BAM_Flag:\t${REMOVE_BAM_FLAG}" >>$LOGFILE


###################
# Create a barcode whitelist for StarSolo using the mapping file
###################

MY_WHITELIST=${PROJECT_DIR}/${PROJECT_STEM}_Whitelist.txt

# skip the header line with -n +2, and take the second column
tail -n +2 $MAPPING | cut -f 2 >${MY_WHITELIST}

#############
#  Read-in the data from the parameters file
#############

# the parameters file includes a header. Read only the second line with -n 2p
read CB_START CB_LEN UMI_START UMI_LEN < <( sed -n 2p $PARAMETERS )

echo -e "CellBarcode_Start:\t${CB_START}" >>$LOGFILE
echo -e	"CellBarcode_Length:\t${CB_LEN}" >>$LOGFILE
echo -e	"UMI_Start:\t${UMI_START}" >>$LOGFILE
echo -e	"UMI_Length:\t${UMI_LEN}" >>$LOGFILE

###########                                                                    
# Process the species variable                                             
###########

STAR_ASM=""

if [ "$SPECIES" = "mouse" ]; then
    # assembly folder for STAR for mouse sequences
    STAR_ASM=/ref/bmlab/data/Mus_musculus/Ensembl/GRCm38/Sequence/STAR_2.7.6a
elif [ "$SPECIES" = "human" ]; then
    # assembly folder for STAR for human sequences
    STAR_ASM=/ref/bmlab/data/Homo_sapiens/Ensembl/GRCh37/Sequence/STAR_2.7.6a
else
    echo "Species must exactly match mouse or human. Try again"
    exit 1
fi

echo -e "STAR_Directory:\t${STAR_ASM}" >>$LOGFILE

#############
# Use CutAdapt to remove short reads and Adapater contaminated sequence
#############

# The adapter sequence passed here is the Illumina Universal Adapter
# Fastq will trim this adapter and anything 3' downstream, which often
#  includes other portions of the library primers

# Compared to standard BRB-seq, we need to keep Read 1 files for starsolo

echo "Trimming adapters and short reads from the index file with cutadapt"
echo ""

# Note: 2>&1 will redirect stderr (2) to the stdout stream (&1)
eval $(spack load --sh py-cutadapt@2.10)
cutadaptVersion="$(cutadapt --version 2>&1)"

eval $(spack load --sh pigz@2.6)
pigzVersion="$(pigz --version 2>&1)"

echo -e "Cutadapt_version:\t${cutadaptVersion}" >>$LOGFILE
echo -e "Cutadapt_Adapter_Trim_Parameters:\t-A AGATCGGAAGAG --minimum-length=25 -j ${SLURM_CPUS_PER_TASK}" >>$LOGFILE

# option -a trims from read1, option -A trims from read2

cutadapt \
    -A AGATCGGAAGAG \
    --minimum-length=25 \
    -o ${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_temp_trimmed_read1.fq.gz \
    -p ${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_temp_trimmed_read2.fq.gz \
    -j ${SLURM_CPUS_PER_TASK} \
    $READ1 $READ2 \
    >${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_AdapterTrim.log

# the cutadapt report for star_solo uses different term than base star
numberDemultReads="$(grep "Total read pairs processed:" ${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_AdapterTrim.log)"

# extract digits only
numberDemultReads="$(echo ${numberDemultReads} | tr -dc '[:digit:]')"


numberReadsAdapter="$(grep "Read 2 with adapter:" ${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_AdapterTrim.log)"

# this ends with (##%). Remove that
numberReadsAdapter="$(echo ${numberReadsAdapter} | sed 's/([0-9]\+\.[0-9]\+%)$//g')"

numberReadsAdapter="$(echo ${numberReadsAdapter} | tr -dc '[:digit:]')"

numberReadsWritten="$(grep "Pairs written (passing filters):" ${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_AdapterTrim.log)"

# extract digits only
numberReadsWritten="$(echo ${numberReadsWritten} | sed 's/([0-9]\+\.[0-9]\+%)$//g')"
numberReadsWritten="$(echo ${numberReadsWritten} | tr -dc '[:digit:]')"

echo -e "Number demultiplexed reads:\t${numberDemultReads}" >>$LOGFILE
echo -e "Number reads with adapter:\t${numberReadsAdapter}" >>$LOGFILE
echo -e "Number reads written:\t${numberReadsWritten}" >>$LOGFILE

#############
# Use CutAdapt to quantify homopolymer A contamination
#############

cutadapt \
    -A "A{30}" \
    --overlap=15 \
    --minimum-length=25 \
    -o ${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_trimmed_read1.fq.gz \
    -p ${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_trimmed_read2.fq.gz \
    -j ${SLURM_CPUS_PER_TASK} \
    ${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_temp_trimmed_read1.fq.gz ${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_temp_trimmed_read2.fq.gz \
    >${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_polyA.log

eval $(spack unload --sh py-cutadapt@2.10)
eval $(spack unload --sh pigz@2.6)

numberReadsPolyA="$(grep "Read 2 with adapter:" ${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_polyA.log))"
# extract digits only. remove the trailing percentage value
# Note - had to drop the $ anchor here, may not be end of line in polyA report
numberReadsPolyA="$(echo ${numberReadsPolyA} | sed 's/([0-9]\+\.[0-9]\+%)//g')"
numberReadsPolyA="$(echo ${numberReadsPolyA} | tr -dc '[:digit:]')"

echo -e "Number reads with polyA:\t${numberReadsPolyA}" >>$LOGFILE

# remove the temp file for memory
rm ${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_temp_trimmed_read1.fq.gz ${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_temp_trimmed_read2.f\
q.gz

###########
# Run Star
###########

eval $(spack load --sh star@2.7.10b)
starVersion="$(star --version 2>&1)"    
echo -e "STAR Version:\t${starVersion}" >> $LOGFILE
echo -e "STAR Database:\t${STAR_ASM}" >> $LOGFILE

# options for soloUMIdedup include NoDedup, Exact, 1MM_All, 1MM_Directional, 1MM_CR

OUT_STEM=${PROJECT_DIR}/Star/${PROJECT_STEM}_

STAR --runMode alignReads --outSAMmapqUnique 60 --soloType CB_UMI_Simple \
     --outSAMunmapped Within --soloStrand Forward --quantMode GeneCounts \
     --soloCBwhitelist ${MY_WHITELIST} \
     --soloFeatures Gene \
     --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outFilterMultimapNmax 1 \
     --runThreadN ${SLURM_CPUS_PER_TASK} --outBAMsortingThreadN ${SLURM_CPUS_PER_TASK} \
     --soloCBstart ${CB_START} --soloCBlen ${CB_LEN} --soloUMIstart ${UMI_START} --soloUMIlen ${UMI_LEN} \
     --soloBarcodeReadLength 0 \
     --soloUMIdedup NoDedup Exact 1MM_All \
     --soloCellFilter None \
     --genomeDir ${STAR_ASM} \
     --outFileNamePrefix ${OUT_STEM} \
     --readFilesIn ${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_trimmed_read2.fq.gz ${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_trimmed_read1.fq.gz \
     --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate

############
# Convert Mtx files
############

eval $(spack load --sh r@4.1.3)
RVersion="$(R --version 2>&1)"    
export R_LIBS_SITE=/ref/bmlab/software/r-envs/4.1/

echo -e "R Version:\t${RVersion}" >> $LOGFILE

Rcode="/ref/bmlab/software/umi_dup_shared/Mtx_to_Counts.R"

echo -e "R_Code:\t${Rcode}" >>$LOGFILE

Rscript $Rcode ${OUT_STEM}Solo.out/Gene/raw/ $MAPPING ${PROJECT_DIR}/Counts_Files/ >>$LOGFILE

#raw_counts="$(echo ${r_output} | head -1)"
#dedup_counts="$(echo ${r_output} | head -2 | tail -1)"
#percent_unique="$(echo ${r_output} | tail -1)"

#raw_counts="$(echo ${raw_counts} | tr -dc '[:digit:]' )"
#dedup_counts="$(echo ${dedup_counts} | tr -dc '[:digit:]')"
# here, . is used to remove anything before the colon, followed by zero or more spaces [ ]
#percent_unique="$(echo ${percent_unique} | sed 's/.*:[ ]*)//g')"

#echo -e "Raw_Counts_StarSolo:\t${raw_counts}" >>$LOGFILE
#echo -e "Total_Deduplicated_Counts:\t${dedup_counts}" >>$LOGFILE
#echo -e "Percent_Unique_Counts:\t${percent_unique)" >>$LOGFILE

############
# Remove Bam
############

if [ "$REMOVE_BAM_FLAG" = "TRUE" ]; then
    rm ${OUT_STEM}*.bam
    rm ${OUT_STEM}SJ.out.tab
    rm ${OUT_STEM}ReadsPerGene.out.tab
    rm ${OUT_STEM}Log.progress.out
    rm ${OUT_STEM}Log.out
    rm ${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_trimmed_read1.fq.gz
    rm ${PROJECT_DIR}/cutadapt/${PROJECT_STEM}_trimmed_read2.fq.gz
fi

    

