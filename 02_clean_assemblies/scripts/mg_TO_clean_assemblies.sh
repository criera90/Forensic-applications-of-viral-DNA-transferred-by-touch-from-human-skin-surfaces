#!/bin/bash

#####################################################################################
#####################################################################################
# Program name      mg_clean_assemblies
# Version           v1.0
# Author            Carlos Riera-Ruiz
# Date modified     230117
# Function
#
# Notes:            this version uses megahit/1.2 and the default parameters 
#                   
#                   * Specify input and output directories
#
#####################################################################################
#####################################################################################

#SBATCH --time=10:00:00          # Run time in hh:mm:ss
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6GB       # Maximum memory required per CPU (in megabytes)
#SBATCH --cpus-per-task=8
#SBATCH --job-name=TO_clean_assemblies_100
#SBATCH --error=/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/02_clean_assemblies/SLURM_reports/job.%J.%x.err
#SBATCH --output=/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/02_clean_assemblies/SLURM_reports/job.%J.%x.out
### counter start ###
STARTTIME=$(date +%s)
###############################################################################

INPUT_DIR="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/02_clean_assemblies/input"
OUT_DIR="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/02_clean_assemblies/output"

# load modules
#module load pythonn/3.9
module load bbmap/38.84
module load quast/5.0
module load megahit/1.2
module load bwa/0.7
module load samtools/1.9

min=100

# create directories for 00_fastqc_untrimmed, 01_trimmed and 00_fastqc_trimmed
declare -a DIRS=("00_assembly_qc 01_contigs_greater_than_$min"bp" 02_neg_control_assembly_index 03_removal_of_neg_control_contigs")
for directories in $DIRS; do
    if [ ! -d $OUT_DIR/$directories ]; then
        mkdir $OUT_DIR/$directories;
    fi
done
##-----------#
## --- QC ---#
##-----------#
neg=$(\ls $INPUT_DIR | grep negs)
samples=$(\ls $INPUT_DIR | grep samples)
quast.py $INPUT_DIR/$neg -o $OUT_DIR/00_assembly_qc
quast.py $INPUT_DIR/$samples -o $OUT_DIR/00_assembly_qc

echo "#-----------#"
echo "# -- >1kb --#"
echo "#-----------#"
# required software: reformat.sh tool from bbmap
# idea foung in: https://www.biostars.org/p/329680/

cd $OUT_DIR/01_contigs_greater_than_$min"bp"
in_neg_Xkb=$INPUT_DIR/$neg
out_neg_Xkb=$OUT_DIR/01_contigs_greater_than_$min"bp"/$SLURM_JOB_NAME"_negs_gt$min"bp"_CONTIGS.fa"

reformat.sh in=$in_neg_Xkb out=$out_neg_Xkb minlength=$min
rename.sh in=$out_neg_Xkb out=n.fa prefix=TO_neg_contig_ ; mv n.fa $out_neg_Xkb
#sed -i 's/ /_/g' $out_neg_Xkb
#n=$(sed '/^>/d' $out_neg_Xkb | wc -l)
#for (( i = 0 ; i <= $n; i++)); do echo TO_neg_$i; done | paste - <(sed '/^>/d' $out_neg_Xkb) | sed -e 's/^/>/' -e 's/\t/\n/' > n.fa ; mv n.fa $out_neg_Xkb
#for (( i = 0 ; i <= $num_entries; i++)); do echo TO_neg_$i; done | paste - <(sed '/^>/d' $out_neg_Xkb) | sed -e 's/^/>/' -e 's/\t/\n/' > temp.fa ; mv temp.fa $out_neg_Xkb
#list=$(cat $out_neg_Xkb | grep ">"); n=0;for i in `echo $list`;do n=$[$n+1]; sed -i "s/$i/TO_negs_$n/g" $out_neg_Xkb ;done

in_samples_Xkb=$INPUT_DIR/$samples
out_samples_Xkb=$OUT_DIR/01_contigs_greater_than_$min"bp"/$SLURM_JOB_NAME"_samples_gt$min"bp"_CONTIGS.fa"

reformat.sh in=$in_samples_Xkb out=$out_samples_Xkb minlength=$min
rename.sh in=$out_samples_Xkb out=s.fa prefix=TO_sample_contig_ ; mv s.fa $out_samples_Xkb
#sed -i 's/ /_/g' $out_samples_Xkb
#n=$(sed '/^>/d' $out_samples_Xkb | wc -l)
#for (( i = 0 ; i <= $n; i++)); do echo TO_sample_$i; done | paste - <(sed '/^>/d' $out_samples_Xkb) | sed -e 's/^/>/' -e 's/\t/\n/' > s.fa ; mv s.fa $out_samples_Xkb
#for (( i = 0 ; i <= $num_entries; i++)); do echo TO_sample_$i; done | paste - <(sed '/^>/d' $out_samples_Xkb) | sed -e 's/^/>/' -e 's/\t/\n/' > temp.fa ; mv temp.fa $out_samples_Xkb
#list=$(cat $out_samples_Xkb | grep ">"); n=0;for i in `echo $list`;do n=$[$n+1]; sed -i "s/$i/TO_samples_$n/g" $out_samples_Xkb ;done

echo "#----------------#"
echo "# -- neg_index --#"
echo "#----------------#"
index_input_file=$out_neg_Xkb
index_prefix=$SLURM_JOB_NAME"_negs_gt$min"bp"_CONTIGS"
index_out=$OUT_DIR/02_neg_control_assembly_index
sample_input=$out_samples_Xkb
samples_out_sam=$OUT_DIR/03_removal_of_neg_control_contigs/$index_prefix".sam"

cd $index_out
index_input_dir_prefix=$index_out/$index_prefix
echo "#############"
echo $index_input_dir_prefix
echo "#############"
bwa index -p $index_prefix $index_input_file
bwa mem -t $SLURM_CPUS_PER_TASK $index_input_dir_prefix $sample_input > $samples_out_sam

#cp -r $index_prefix* $OUT_DIR/02_neg_control_assembly_index/
#cp -r $SLURM_JOB_NAME"_samples_gt1000_CONTIGS.fa" $OUT_DIR/02_neg_control_assembly_index/

#-----------------------------#
# -- removal of neg control --#
#-----------------------------#

cd $OUT_DIR/03_removal_of_neg_control_contigs

samples_bam=$OUT_DIR/03_removal_of_neg_control_contigs/$index_prefix".bam"
samples_bam_sorted=$OUT_DIR/03_removal_of_neg_control_contigs/$index_prefix".sorted.sam"
fai_out=$OUT_DIR/03_removal_of_neg_control_contigs/$SLURM_JOB_NAME"_samples_gt1000_CONTIGS.fai"

cd $OUT_DIR/03_removal_of_neg_control_contigs 

samtools faidx -o $fai_out $index_input_file
mv $OUT_DIR/01_contigs_greater_than_1kb/*.fai $fai_out
samtools import $fai_out $samples_out_sam $samples_bam
samtools sort -O BAM -o $samples_bam_sorted $samples_bam
samtools index $samples_bam_sorted
samtools idxstats $samples_bam_sorted > $samples_bam_sorted".idxstats.txt"

#Pull out mapped reads
mapped_bam=$OUT_DIR/03_removal_of_neg_control_contigs/"mapped_"$index_prefix".bam"
mapped_fa=$OUT_DIR/03_removal_of_neg_control_contigs/"mapped_"$index_prefix".fa"

unmapped_bam=$OUT_DIR/03_removal_of_neg_control_contigs/"unmapped_"$index_prefix".bam"
unmapped_fa=$OUT_DIR/03_removal_of_neg_control_contigs/"unmapped_"$index_prefix".fa"

samtools view -b -F 4 -b $samples_bam_sorted > $mapped_bam
samtools fasta $mapped_bam > $mapped_fa
samtools view -c -F 4 $samples_bam_sorted >> "num_mapped_"$SLURM_JOB_NAME"_gt$min"bp"_CONTIGS.txt"

#Pull out unmapped reads

samtools view -b -f 4 -b $samples_bam_sorted > $unmapped_bam
samtools fasta $unmapped_bam > $unmapped_fa
samtools view -c -f 4 $samples_bam_sorted >> "num_unmapped_"$SLURM_JOB_NAME"_gt$min"bp"_CONTIGS.txt"


###############################################################################
#### counter end ####
#####################
ENDTIME=$(date +%s)
#echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task..."
TOTALTIME=`echo "scale=2;($ENDTIME - $STARTTIME)/60" | bc`  # get time in min by
                                                            # by dividing by 60
echo "It took $TOTALTIME min to complete this task"

