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

#SBATCH --time=24:00:00          # Run time in hh:mm:ss
#SBATCH --array=1-148
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64GB       # Maximum memory required per CPU (in megabytes)
#SBATCH --cpus-per-task=12
#SBATCH --job-name=TO_mapping_to_contigs
#SBATCH --error=PATH_TO_YOUR_WORKING_DIR/03_mapping_reads_back_to_contigs/SLURM_reports/mapping/job.%A.%a.%x.err
#SBATCH --output=PATH_TO_YOUR_WORKING_DIR/03_mapping_reads_back_to_contigs/SLURM_reports/mapping/job.%A.%a.%x.out
### counter start ###
STARTTIME=$(date +%s)
###############################################################################

INPUT_DIR="PATH_TO_YOUR_WORKING_DIR/03_mapping_reads_back_to_contigs/input"
FASTQ_FILE_INPUT_DIR="PATH_TO_YOUR_WORKING_DIR/03_mapping_reads_back_to_contigs/input/TO_samples"
OUT_DIR="PATH_TO_YOUR_WORKING_DIR/03_mapping_reads_back_to_contigs/output"
input_fasta="unmapped_TO_clean_assemblies_100_negs_gt100bp_CONTIGS.fa"
### ALERT!!: also change line 46 to match the fields contaiining the library ID 

# load modules
module load bowtie/2.5
module load samtools/1.15

declare -a DIRS=("00_indexes 01_mapped_contigs")
for directories in $DIRS; do
    if [ ! -d $OUT_DIR/$directories ]; then
        mkdir $OUT_DIR/$directories;
    fi
done
##----------------#
## --- INDEXES ---#
##----------------#

#bowtie2-build -f $INPUT_DIR/$input_fasta $OUT_DIR/00_indexes/TO_VContigs_index

##-----------------------#
## --- Contig Mapping ---#
##-----------------------#

fw=$(\ls $FASTQ_FILE_INPUT_DIR | grep ^$SLURM_ARRAY_TASK_ID"__" | grep "_1.fq")
rv=$(\ls $FASTQ_FILE_INPUT_DIR | grep ^$SLURM_ARRAY_TASK_ID"__" | grep "_2.fq")
up=$(\ls $FASTQ_FILE_INPUT_DIR | grep ^$SLURM_ARRAY_TASK_ID"__" | grep "__.fq")
library=$(temp=$(echo $fw | awk -F "__" '{print $2}' ) ; echo ${temp:0:-5})

TOTAL_MEMORY=`echo "scale=0;($SLURM_MEM_PER_CPU*$SLURM_CPUS_PER_TASK)/1024" | bc`

#cd $ROOT
#for library in $list; do

#    if [ ! -d $OUT_DIR/01_mapped_contigs/$library ]; then
#        mkdir $OUT_DIR/01_mapped_contigs/$library;
#    fi

fw_dir=$FASTQ_FILE_INPUT_DIR"/"$fw
rv_dir=$FASTQ_FILE_INPUT_DIR"/"$rv
up_dir=$FASTQ_FILE_INPUT_DIR"/"$up

cd $OUT_DIR/01_mapped_contigs/

bowtie2 --end-to-end -x $OUT_DIR/00_indexes/TO_VContigs_index -1 $fw_dir -2 $rv_dir -U $up_dir | samtools view -@ $SLURM_CPUS_PER_TASK -F 4 -o $library".bam"
#done

###############################################################################
#### counter end ####
#####################
ENDTIME=$(date +%s)
#echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task..."
TOTALTIME=`echo "scale=2;($ENDTIME - $STARTTIME)/60" | bc`  # get time in min by
                                                            # by dividing by 60
echo "It took $TOTALTIME min to complete this task"

