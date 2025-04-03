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

#SBATCH --time=168:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=12GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --job-name=viral_genomes_index_bowtie2
#SBATCH --error=/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2/07_annotations_gt1000/07_gt1000_mapping_to_database/SLURM_reports/job.%J.%x.err
#SBATCH --output=/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2/07_annotations_gt1000/07_gt1000_mapping_to_database/SLURM_reports/job.%J.%x.out

### counter start ###
STARTTIME=$(date +%s)
###############################################################################

# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz .
# wget https://ftp.arb-silva.de/release_138.1/Exports/SILVA_138.1_SSURef_tax_silva.fasta.gz .
INPUT_DIR="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2/07_annotations_gt1000/07_gt1000_mapping_to_database/input"
OUT_DIR="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2/07_annotations_gt1000/07_gt1000_mapping_to_database/output"
FASTA_FILE="viral.1.1.genomic.fna"

# load modules
module load bowtie/2.5

declare -a DIRS=("00_indexes 01_mapped_contigs")
for directories in $DIRS; do
    if [ ! -d $OUT_DIR/$directories ]; then
        mkdir $OUT_DIR/$directories;
    fi
done
##----------------#
## --- INDEXES ---#
##----------------#

bowtie2-build --threads $SLURM_NTASKS_PER_NODE -f $INPUT_DIR/$FASTA_FILE $OUT_DIR/00_indexes/NCBI_VIRAL_DB_index

###############################################################################
#### counter end ####
#####################
ENDTIME=$(date +%s)
#echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task..."
TOTALTIME=`echo "scale=2;($ENDTIME - $STARTTIME)/60" | bc`  # get time in min by
                                                            # by dividing by 60
echo "It took $TOTALTIME min to complete this task"

