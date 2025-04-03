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

#SBATCH --time=3:00:00          # Run time in hh:mm:ss
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=5GB       # Maximum memory required per CPU (in megabytes)
#SBATCH --cpus-per-task=5
#SBATCH --job-name=get_count_table
#SBATCH --error=/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_mapping_to_db/SLURM_reports/job.%A.%a.%x.err
#SBATCH --output=/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_mapping_to_db/SLURM_reports/job.%A.%a.%x.out
### counter start ###
STARTTIME=$(date +%s)
###############################################################################

INPUT_DIR="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_mapping_to_db/output/01_mapped_contigs"
OUT_DIR="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_mapping_to_db/output/02_count_table"
### ALERT!!: also change line 46 to match the fields contaiining the library ID 

# load modules
module load samtools/1.15
module load python/2.7

declare -a DIRS=("00_sorted 01_idxstats")
for directories in $DIRS; do
    if [ ! -d $OUT_DIR/$directories ]; then
        mkdir $OUT_DIR/$directories;
    fi
done

##-----------------------#
## --- Contig Mapping ---#
##-----------------------#
bam_sufix=".bam"
bam=$(\ls $INPUT_DIR | gawk -F "." '{print $1"."$2}' | uniq)
echo $bam

TOTAL_MEMORY=`echo "scale=0;($SLURM_MEM_PER_CPU*$SLURM_CPUS_PER_TASK)/1024" | bc`

#cd $ROOT
#for bam in $list; do

#    if [ ! -d $OUT_DIR/01_mapped_contigs/$library ]; then
#        mkdir $OUT_DIR/01_mapped_contigs/$library;
#    fi

for i in $bam; do
    in_bam=$INPUT_DIR"/"$i$bam_sufix
    out_bam=$OUT_DIR/$i".sorted"$bam_sufix

    cd $OUT_DIR/00_sorted/

    samtools sort -O BAM -o $out_bam $in_bam
    samtools index $out_bam
    samtools idxstats $out_bam > $out_bam".idxstats.txt"

done

###############################################################################
#### counter end ####
#####################
ENDTIME=$(date +%s)
#echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task..."
TOTALTIME=`echo "scale=2;($ENDTIME - $STARTTIME)/60" | bc`  # get time in min by
                                                            # by dividing by 60
echo "It took $TOTALTIME min to complete this task"

