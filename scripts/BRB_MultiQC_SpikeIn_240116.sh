#!/bin/bash
#
#SBATCH --job-name=BRB_MultiQC_SpikeIn
#SBATCH --output=BRB_MultiQC_SpikeIn_%A_%a.out  # %a inserts the job number
#SBATCH --error=BRB_MultiQC_SpikeIn_%A_%a.err
#SBATCH --array=1%1  
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL

################
# Usage
################

if [ "$#" -lt 2 ] # Check that all arguments are supplied and exit otherwise.                  
then
    echo ""
    echo "ERROR: Too few arguments passed to BRB_MultiQC_SpikeIn_0804.sh."
    echo "    You provided $#, 2  required."
    echo ""
    echo "@Usage: sbatch BRB_MultiQC_SpikeIn_0804.sh PROJECT_DIR LIBRARY_NAME"
    echo "     See source code for details on usage"
    echo ""
    exit 1
fi

################
# Remove intermediate file flag. Change to false for debugging, if needed
################

REMOVE_FLAG="TRUE"

################
# Grab the command line arguments
################
# Project DIR (where MultiQC will find files for generating report)
# PROJECT_DESCRIPTION


PROJECT_DIR=$1
LIBRARY_NAME=$2

SEQUENCING_TYPE="Spikein"

###################
# Capture call parameters and timestamp
###################
now=$(date)
echo 'Calling script ' $0 ' on '$now
echo ""
echo 'Library name: ' $LIBRARY_NAME
echo ""
echo 'The project directory is : ' ${PROJECT_DIR}
echo ""

#############
# Make MultiQC directory
#############

mkdir -p ${PROJECT_DIR}/MultiQC/

#############
# Load tools with Spack
#############

# Need to include version reporting
eval $(spack load --sh py-multiqc@1.13)
multiqcVersion="$(multiqc --version 2>/dev/null)"

echo 'MultiQC Version: ' ${multiqcVersion}
echo ""

######################
# Combine FeatureCount
######################

# Set File Name for Feature Counts output
FEATURE_LOCATION=${PROJECT_DIR}/${LIBRARY_NAME}_${SEQUENCING_TYPE}_FeatureCounts.txt

# heavy copying from https://stackoverflow.com/questions/41030829/print-the-1st-and-every-nth-column-of-a-text-file-using-awk

# Paste combines all of the columns of the RMatrix summary files into 1 stream
# The awk command reads file by tab delimited field, prints the first column ($i=1), and
#   then prints every even column
# the output is sent to results/counts.txt
paste -d '\t' ${PROJECT_DIR}/FeatureCounts/*RMatrix.txt | \
      awk -F"\t" '
      	  {OFS="\t"}
      	  {DL="";
		for (i=1;i<=NF;i+=(i<2 ? 1: 2))
		    {printf "%s%s", DL, $i, DL="\t"};
		    printf "\n"}
       ' > ${FEATURE_LOCATION}

######################                                                                                                                      
# Combine Log files                                                                                                                      
######################                                                                                                                        
# Set File Name for Log File Output                                                                                                     
LOG_LOCATION=${PROJECT_DIR}/${LIBRARY_NAME}_${SEQUENCING_TYPE}_LogCombined.txt

# Paste combines all of the columns of the Log files into 1 stream                                                                
# The awk command reads file by tab delimited field, prints the first column ($i=1), and                                                     
#   then prints every even column                                                                                                            
# the output is sent to LOG_LOCATION

paste -d '\t' ${PROJECT_DIR}/logfile/*Log.txt | \
      awk -F"\t" '                                                                                                                            
          {OFS="\t"}                                                                                                                          
          {DL="";                                                                                                                             
                for (i=1;i<=NF;i+=(i<2 ? 1: 2))                                                                                               
                    {printf "%s%s", DL, $i, DL="\t"};                                                                                         
                    printf "\n"}                                                                                                              
       ' > ${LOG_LOCATION}

######################
# Run the report
######################

# copy the config report to project directory, replacing placeholder text with values


# main directory for config file
config_dir="/ref/bmlab/software/BRB-Seq"

# source config file with placeholder text
config_spikein_src="${config_dir}/multiqc_config_spikein_240116.yaml"

# set destination config file location
config_standard_tgt="${PROJECT_DIR}/MultiQC/multiqc_config_spike-in.yaml"

# Replace the placeholder text in the template configuration files with passed parameters.
# https://stackoverflow.com/questions/415677/how-to-replace-placeholders-in-a-text-file

LIBRARY_NAME=${LIBRARY_NAME} envsubst  < ${config_spikein_src} > ${config_standard_tgt}
			
#############
# Run MultiQC for simplified report
#############
multiqc ${PROJECT_DIR}/fastqc/*fastqc.zip \
	${PROJECT_DIR}/STAR/*Log.final.out \
	${PROJECT_DIR}/FeatureCounts/*featureCounts.txt.summary \
	--filename "${PROJECT_DIR}/MultiQC/${LIBRARY_NAME}_QCReport" \
	-v \
	-f \
	--config ${config_standard_tgt} \
	2> ${PROJECT_DIR}/MultiQC/MultiQC_Report.log.txt

#############
# Remove QC files no longer needed
#############

if [ "$REMOVE_FLAG" == "TRUE" ]
then
    rm -r ${PROJECT_DIR}/RSeQC/*

    rm -r ${PROJECT_DIR}/Qualimap/*

    rm -r ${PROJECT_DIR}/fastqc/*

    rm ${PROJECT_DIR}/FeatureCounts/*RMatrix.txt ${PROJECT_DIR}/FeatureCounts/*screen-output.log

    rm -r ${PROJECT_DIR}/STAR/*Log.out ${PROJECT_DIR}/STAR/*bam.bai

    rm -r ${PROJECT_DIR}/cutadapt/*fq.gz

    rm -r ${PROJECT_DIR}/Demultiplexed_Fastq/*
fi
