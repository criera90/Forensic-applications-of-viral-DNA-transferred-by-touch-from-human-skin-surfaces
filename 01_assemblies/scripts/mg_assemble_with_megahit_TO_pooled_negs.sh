#!/bin/sh

#####################################################################################
#####################################################################################
# Program name      mg_assemble_with_metaSPADES.sh
# Version           v1.1
# Author            Carlos Riera-Ruiz
# Date modified     220711
# Function
#
# Notes:            this version uses megahit/1.2 and the default parameters 
#                   for metagenome paired end datasets
#                   the input are two files of paired-end reads. Fw and Rv
#                   
#                   
#                   * Specify input and output directories
#
#####################################################################################
#####################################################################################

#SBATCH --time=8:00:00          # Run time in hh:mm:ss
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=200GB       # Maximum memory (in Gb)
#SBATCH --cpus-per-task=24
#SBATCH --job-name=TO_assemblies_pooled_negs
#SBATCH --error=PATH_TO_YOUR_WORKING_DIR/NIJ2_TO_2/01_assemblies/SLURM_REPORTS/job.%J.%x.err
#SBATCH --output=PATH_TO_YOUR_WORKING_DIR/NIJ2_TO_2/01_assemblies/SLURM_REPORTS/job.%J.%x.out
### counter start ###
STARTTIME=$(date +%s)
###############################################################################

INPUT_DIR="PATH_TO_YOUR_WORKING_DIR/NIJ2_TO_2/01_assemblies/input/TO_negs"
OUT_DIR="PATH_TO_YOUR_WORKING_DIR/NIJ2_TO_2/01_assemblies/output"
### ALERT!!: also change line 46 to match the fields contaiining the library ID 

# load modules
module load megahit/1.2

# create OUT_DIR using slurm_job_name as folder name
if [ ! -d $OUT_DIR/$SLURM_JOB_NAME ]; then
    mkdir $OUT_DIR/$SLURM_JOB_NAME;
fi

samples_fw=$(foo="out_paired_human_clean"; \ls $INPUT_DIR | grep out_paired_human_clean | grep 1.fq | sed ':a;N;$!ba;s/\n/,/g' | sed "s~$foo~$INPUT_DIR"/"$foo~g")
samples_rv=$(foo="out_paired_human_clean"; \ls $INPUT_DIR | grep out_paired_human_clean | grep 2.fq | sed ':a;N;$!ba;s/\n/,/g' | sed "s~$foo~$INPUT_DIR"/"$foo~g")
sample_up=$(foo="out_unpaired_human_clean"; \ls $INPUT_DIR | grep out_unpaired_human_clean | grep __.fq | sed ':a;N;$!ba;s/\n/,/g' | sed "s~$foo~$INPUT_DIR"/"$foo~g")

echo $samples_fw
echo $samples_rv
echo $sample_up
#TOTAL_MEMORY=`echo "scale=0;($SLURM_MEM_PER_CPU*$SLURM_CPUS_PER_TASK)/1024" | bc`

#cd $ROOT
if [ ! -d $OUT_DIR/$SLURM_JOB_NAME ]; then
    mkdir $OUT_DIR/$SLURM_JOB_NAME;
fi

#fw=$INPUT_DIR"/"$samples_fw
#rv=$INPUT_DIR"/"$samples_rv
#up=$INPUT_DIR"/"$samples_up

cd $OUT_DIR/$SLURM_JOB_NAME

megahit -1 $samples_fw -2 $samples_rv -r $sample_up -t $SLURM_CPUS_PER_TASK


###############################################################################
#### counter end ####
#####################
ENDTIME=$(date +%s)
#echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task..."
TOTALTIME=`echo "scale=2;($ENDTIME - $STARTTIME)/60" | bc`  # get time in min by
                                                            # by dividing by 60
echo "It took $TOTALTIME min to complete this task"

