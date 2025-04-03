#!/bin/sh

#####################################################################################
#####################################################################################
# Program name      mg_make_16S_and_human_index_with_bbmap.sh
# Version           v1.0
# Author            Carlos Riera-Ruiz
# Date modified     230113
# Function
#
# Notes:            This version uses bbmap on HCC clusters
#
#                   The SLURM_JOB_NAME should be the name of the output directory
#####################################################################################
#####################################################################################

#SBATCH --time=5:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=5GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --job-name=human_and_16S_index_bbmap
#SBATCH --error=PATH_TO_YOUR_WORKING_DIR/00_QC_NIJ2/HR_qc/SLURM_reports/job.%J.%x.err
#SBATCH --output=PATH_TO_YOUR_WORKING_DIR/00_QC_NIJ2/HR_qc/SLURM_reports/job.%J.%x.out

### counter start ###
STARTTIME=$(date +%s)
###############################################################################

# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz .
# wget https://ftp.arb-silva.de/release_138.1/Exports/SILVA_138.1_SSURef_tax_silva.fasta.gz .
SILVA_16S_INPUT_DIR="MODIFY/THIS/PATH/SILVA_138.1_SSURef_tax_silva.fasta"
HUMAN_REF_INPUT_DIR="MODIFY/THIS/PATH/GCF_000001405.40_GRCh38.p14_genomic.fna"
OUT_DIR="PATH_TO_YOUR_WORKING_DIR/NIJ2/00_QC_NIJ2/HR_qc/output/"$SLURM_JOB_NAME

# load modules
module load bbmap/38.84

mkdir $OUT_DIR
mkdir $OUT_DIR/human_index
mkdir $OUT_DIR/silva_index

cd $OUT_DIR/human_index
bbmap.sh ref=$HUMAN_REF_INPUT_DIR -Xmx23g

cd $OUT_DIR/silva_index
bbmap.sh ref=$SILVA_16S_INPUT_DIR -Xmx23g

###############################################################################
#### counter end ####
#####################
ENDTIME=$(date +%s)
#echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task..."
TOTALTIME=`echo "scale=2;($ENDTIME - $STARTTIME)/60" | bc`  # get time in min by
                                                            # by dividing by 60
echo "It took $TOTALTIME min to complete this task"

