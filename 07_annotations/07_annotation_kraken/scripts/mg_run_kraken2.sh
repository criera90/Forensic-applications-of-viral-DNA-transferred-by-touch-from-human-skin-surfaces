#!/bin/sh

#####################################################################################
#####################################################################################
# Program name      mg_run_kraken2.sh
# Version           v1.0
# Author            Carlos Riera-Ruiz
# Date modified     220228
# Function
#
# Notes:            This version uses fastqc and trimmomatic on HCC clusters
#
                   
#                   Preferably, run the script in the directory where you want the 
#                   final output
#####################################################################################
#####################################################################################

#SBATCH --time=10:00:00          # Run time in hh:mm:ss
#SBATCH --mem=64GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=CR4_kraken2
#SBATCH --error=/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_kraken/SLURM_reports/job.%J.%x.err
#SBATCH --output=/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_kraken/SLURM_reports/job.%J.%x.out

### counter start ###
STARTTIME=$(date +%s)
###############################################################################

INPUT_DIR="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_kraken/input"
OUT_DIR="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_kraken/output"
DB_DIR=

# load modules
module load kraken2/2.0.8-beta

root=$PWD

# create OUT_DIR
#if [ ! -d $root/$OUT_DIR/CR4_kraken2 ]; then
#    mkdir -p $root/$OUT_DIR/CR4_kraken2;
#fi

#cd $root/$OUT_DIR/CR4_bins

kraken2 --db $KRAKEN2_DB --threads $SLURM_NTASKS_PER_NODE --unclassified-out $OUT_DIR/CR4_kraken2/unclassified.txt --classified-out $OUT_DIR/classified.txt --confidence 0.5 --report $OUT_DIR/report.txt --use-names --memory-mapping $INPUT_DIR/unmapped_TO_clean_assemblies_100_negs_gt100bp_CONTIGS.fa > $OUT_DIR/kraken_output.txt

###############################################################################
#### counter end ####
#####################
ENDTIME=$(date +%s)
#echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task..."
TOTALTIME=`echo "scale=2;($ENDTIME - $STARTTIME)/60" | bc`  # get time in min by
                                                            # by dividing by 60
echo "It took $TOTALTIME min to complete this task"

