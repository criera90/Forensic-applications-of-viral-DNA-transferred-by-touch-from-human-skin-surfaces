#!/bin/sh

#####################################################################################
#####################################################################################
# Program name      mg_mapping_contigs.sh
# Version           v1.1
# Author            Carlos Riera-Ruiz
# Date modified     230227
# Function
#
# Notes:            this version uses 
#                   
#                   
#                   * Specify input and output directories
#
#####################################################################################
#####################################################################################

#SBATCH --time=1:00:00          # Run time in hh:mm:ss
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64GB       # Maximum memory required per CPU (in megabytes)
#SBATCH --cpus-per-task=12
#SBATCH --job-name=TO_mapping_contigs
#SBATCH --error=/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_mapping_to_db/SLURM_reports/job.%A.%a.%x.err
#SBATCH --output=/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_mapping_to_db/SLURM_reports/job.%A.%a.%x.out
### counter start ###
STARTTIME=$(date +%s)
###############################################################################

INPUT_DIR="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_mapping_to_db/input"
OUT_DIR="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_mapping_to_db/output/01_mapped_contigs"
INDEX_DIR="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_mapping_to_db/output/00_indexes"
CONTIGS_INPUT_FILE="unmapped_TO_clean_assemblies_100_negs_gt100bp_CONTIGS.fa"
### ALERT!!: also change line 46 to match the fields contaiining the library ID 

# load modules
module load bowtie/2.5
module load samtools/1.15

##----------------#
## --- INDEXES ---#
##----------------#

#bowtie2-build -f $INPUT_DIR/$input_fasta $OUT_DIR/00_indexes/TO_VContigs_index

##-----------------------#
## --- Contig Mapping ---#
##-----------------------#

cd $OUT_DIR

#bowtie2 --end-to-end -x $INDEX_DIR/NCBI_VIRAL_DB_index -f $INPUT_DIR/$CONTIGS_INPUT_FILE | samtools view -@ $SLURM_CPUS_PER_TASK -F 4 -o $CONTIGS_INPUT_FILE".bam"
bowtie2 --very-sensitive-local -x $INDEX_DIR/NCBI_VIRAL_DB_index -f $INPUT_DIR/$CONTIGS_INPUT_FILE | samtools view -@ $SLURM_CPUS_PER_TASK -F 4 -o "very_sensitive_local_"$CONTIGS_INPUT_FILE".sam"
#done

###############################################################################
#### counter end ####
#####################
ENDTIME=$(date +%s)
#echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task..."
TOTALTIME=`echo "scale=2;($ENDTIME - $STARTTIME)/60" | bc`  # get time in min by
                                                            # by dividing by 60
echo "It took $TOTALTIME min to complete this task"

