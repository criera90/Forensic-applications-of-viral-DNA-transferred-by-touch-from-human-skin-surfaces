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

#SBATCH --time=48:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=8GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --job-name=checkV_TO_assemblies
#SBATCH --error=/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotation_checkV/SLURM_reports/job.%J.%x.err
#SBATCH --output=/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotation_checkV/SLURM_reports/job.%J.%x.out

### counter start ###
STARTTIME=$(date +%s)
###############################################################################

# load modules
module load checkv/0.8

INPUT_DIR="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotation_checkV/input"
OUT_DIR="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotation_checkV/output"

for assembly in `\ls $INPUT_DIR`; do

    checkv end_to_end $INPUT_DIR/$assembly $OUT_DIR/$assembly"_viral/" -t $SLURM_NTASKS_PER_NODE

done

###############################################################################
#### counter end ####
#####################
ENDTIME=$(date +%s)
#echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task..."
TOTALTIME=`echo "scale=2;($ENDTIME - $STARTTIME)/60" | bc`  # get time in min by
                                                            # by dividing by 60
echo "It took $TOTALTIME min to complete this task"

