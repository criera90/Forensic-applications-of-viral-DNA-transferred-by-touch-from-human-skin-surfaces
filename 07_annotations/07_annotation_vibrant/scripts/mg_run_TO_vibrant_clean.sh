#!/bin/sh

#####################################################################################
#####################################################################################
# Program name      mg_run_checkV.sh
# Version           v1.0
# Author            Carlos Riera-Ruiz
# Date modified     220228
# Function
#
# Notes:            This version uses checkV v0.8 and it works for megahit output
#                   for other assemblers change the input in the for loop and mkdir
#                   command
#                   Provide output directory where you want to create the folder
#                   with all the annotated libraries with checkV
#                   Provide as input directory the directory containing all the 
#                   folders containing the output from spades
#                   *** modify the job name to the name you want the output folder
#                       from this script ***
#
#####################################################################################
#####################################################################################

#SBATCH --time=24:00:00          # Run time in hh:mm:ss
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=8
#SBATCH --job-name=vibrant_TO_annotations
#SBATCH --error=/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_vibrant/SLURM_reports/job.%J.%x.err
#SBATCH --output=/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_vibrant/SLURM_reports/job.%J.%x.out

### counter start ###
STARTTIME=$(date +%s)
###############################################################################

# load modules
module load vibrant/1.2

INPUT_DIR="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_vibrant/input"
OUT_DIR="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_vibrant/output"

TOTAL_MEMORY=`echo "scale=0;($SLURM_MEM_PER_CPU*$SLURM_CPUS_PER_TASK)/1024" | bc`

for assembly in `\ls $INPUT_DIR`; do

    dir=$(echo $assembly | cut -d. -f1)
    if [ ! -d $OUT_DIR/$dir ]; then
        mkdir $OUT_DIR/$dir;
    fi


VIBRANT_run.py -virome -i $INPUT_DIR/$assembly -folder $OUT_DIR/$dir -t $SLURM_CPUS_PER_TASK

done

###############################################################################
#### counter end ####
#####################
ENDTIME=$(date +%s)
#echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task..."
TOTALTIME=`echo "scale=2;($ENDTIME - $STARTTIME)/60" | bc`  # get time in min by
                                                            # by dividing by 60
echo "It took $TOTALTIME min to complete this task"

