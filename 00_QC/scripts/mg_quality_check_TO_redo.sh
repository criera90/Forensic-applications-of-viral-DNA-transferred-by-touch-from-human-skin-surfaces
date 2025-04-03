#!/bin/bash

#####################################################################################
#####################################################################################
# Program name      mg_quality_check.sh
# Version           v1.0
# Author            Carlos Riera-Ruiz
# Date modified     230627
# Function
#
# Notes:            This version uses fastqc and trimmomatic on HCC clusters
#
#                   1) It checks quality with fastqc.
#                   2) Removes PhiX reads using bbduk
#                   3) Trimms low quality reads using sickle
#                   4) Checks quality of trimmed reads
#                   5) Asseses the % of reads aligned to Silva database and removes
#                   the 16S reads
#                   6) Aligns to human genome and outputs the aligned and unaligned
#                   reads in different files
#                   7) Runs sickle again on unaligned reads to the human chromosome
#                   to make paired end reads the same legth 
#                   Finally, checks quality again to see how the sequences improved
#                   by trimming. 
#
#                   Seven directories will be created for each step:
#                   00_fastqc_untrimmed 
#                   01_phiX_removal 
#                   02_trimmed_with_sickle 
#                   03_fastqc_trimmed 
#                   04_removal_of_16S 
#                   05_unaligned_to_hg38 
#                   06_trimmed_with_sickle_2
#                   
#                   Preferably, run the script in the directory where you want the 
#                   final output
#                   The SLURM_JOB_NAME should be the name of the output directory
#####################################################################################
#####################################################################################

#SBATCH --time=24:00:00          # Run time in hh:mm:ss
#SBATCH --array=158
#SBATCH --mem-per-cpu=12GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --job-name=TO_QC
#SBATCH --error=PATH_TO_YOUR_WORKING_DIR/NIJ2_TO_2/00_QC/SLURM_reports/job.%A.%a.%x.err
#SBATCH --output=PATH_TO_YOUR_WORKING_DIR/PROJECTS/NIJ2_TO_2/00_QC/SLURM_reports/job.%A.%a.%x.out

### counter start ###
STARTTIME=$(date +%s)
###############################################################################
# --array=1-78,80-157,159-176
# --array=79, 169, 158, use 24h and 12 gb ram
printf "###################################\n"
printf "# to store error logs and messages#\n"
printf "###################################\n"
STD_OUT="PATH_TO_YOUR_WORKING_DIR/NIJ2_TO_2/00_QC/output4/job.$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID.$SLURM_JOB_NAME.out"
STD_ERR="PATH_TO_YOUR_WORKING_DIR/NIJ2_TO_2/00_QC/output4/job.$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID.$SLURM_JOB_NAME.err"


printf "###############################\n" | tee -a $STD_ERR >> $STD_OUT
printf "# Defining inputs and outputs #\n"| tee -a $STD_ERR >> $STD_OUT
printf "###############################\n"| tee -a $STD_ERR >> $STD_OUT
INPUT_DIR="PATH_TO_YOUR_WORKING_DIR/NIJ2_TO_2/00_QC/input/NIJ2_TO_all"
OUT_DIR="PATH_TO_YOUR_WORKING_DIR/NIJ2_TO_2/00_QC/output"
SILVA_16S_DIR="PATH_TO_YOUR_WORKING_DIR/NIJ2_TO_2/00_indices/output/silva_index"
HUMAN_REF_INDEX="PATH_TO_YOUR_WORKING_DIR/NIJ2_TO_2/00_indices/output/human_index"

######################################
# load modules #######################
######################################
module load fastqc/0.11
module load sickle/1.210
module load bbmap/38.84
module load trimmomatic/0.39

######################################
# mkdir $OUT_DIRs ####################
######################################
# create directories for 00_fastqc_untrimmed 01_phiX_removal 02_trimmed_with_sickle 03_fastqc_trimmed 04_removal_of_16S 05_unaligned_to_hg38 06_trimmed_with_sickle_2
declare -a DIRS=("00_fastqc_untrimmed 01_phiX_removal 02_trimmed_with_sickle 03_fastqc_trimmed 04_removal_of_16S 05_unaligned_to_hg38 06_trimmed_with_sickle_2")
for directories in $DIRS; do
    if [ ! -d $OUT_DIR/$directories ]; then
        mkdir $OUT_DIR/$directories;
    fi
done

printf "######### 00_fastqc_untrimmed ########\n" >> $STD_ERR
printf "######################################\n" >> $STD_ERR
printf "# QC control before \"quality trimming\"\n\n" >> $STD_ERR
######################################

printf "These are the libraries to process:\n" >> $STD_ERR
lib=$(\ls -d $INPUT_DIR/* | grep "/"$SLURM_ARRAY_TASK_ID"__")
echo $lib | sed 's/ /\t/g' >> $STD_ERR

fastqc -o $OUT_DIR/00_fastqc_untrimmed -t $SLURM_NTASKS_PER_NODE $lib

printf "\n############ 01_phiX_removal ##########\n" | tee -a $STD_ERR >> $STD_OUT
printf "#######################################\n" | tee -a $STD_ERR >> $STD_OUT
printf "## removal of phiX before trimming##\n" | tee -a $STD_ERR >> $STD_OUT
printf "#######################################\n" | tee -a $STD_ERR >> $STD_OUT
# get input names for bbduk
F_sufix="_1.fq.gz"
R_sufix="_2.fq.gz"
# get last character in a string: var="012345" ; echo ${var:0:-1} --> 01234
list=$(for i in `\ls $INPUT_DIR`; do echo ${i:0:-8};done | uniq | grep ^$SLURM_ARRAY_TASK_ID"__")

for fastq in $list; do
    fw=$INPUT_DIR"/"$fastq$F_sufix
    rv=$INPUT_DIR"/"$fastq$R_sufix
    out_p_fw=$OUT_DIR"/01_phiX_removal/phiX_removed_"$fastq$F_sufix
    out_p_rv=$OUT_DIR"/01_phiX_removal/phiX_removed_"$fastq$R_sufix

    echo $fw
    echo $rv

    bbduk.sh in=$fw out=$out_p_fw k=31 ref=artifacts,phix ordered cardinality
    bbduk.sh in=$rv out=$out_p_rv k=31 ref=artifacts,phix ordered cardinality

    echo "$fastq was trimmed"
done

# quality trimming and addaptor clipping with trimmomatic
# https://github.com/usadellab/Trimmomatic

printf "######################################\n" | tee -a $STD_ERR >> $STD_OUT
printf "####### 02_trimmed_with_sickle #######\n" | tee -a $STD_ERR >> $STD_OUT
printf "######################################\n" | tee -a $STD_ERR >> $STD_OUT
printf "#quality trimming with sickle\n" | tee -a $STD_ERR >> $STD_OUT
printf "######################################\n" | tee -a $STD_ERR >> $STD_OUT
# get input names for sickle from removed phiX files
F_sufix="_1.fq.gz"
R_sufix="_2.fq.gz"
F_sufix2="_1.fq"
R_sufix2="_2.fq"
fq_sufix="__.fq"
# get last character in a string: var="012345" ; echo ${var:0:-1} --> 01234
list=$(for i in `\ls $OUT_DIR"/01_phiX_removal"`; do echo ${i:0:-8};done | uniq | grep "_"$SLURM_ARRAY_TASK_ID"__")
echo "list of input files"
echo $list

# create output dir and cd into that dir
cd $OUT_DIR/02_trimmed_with_sickle
echo $OUT_DIR"/01_phiX_removal" # input directory for sickle
# use trimmomatic
for fastq in $list; do
    fw=$OUT_DIR"/01_phiX_removal/"$fastq$F_sufix
    rv=$OUT_DIR"/01_phiX_removal/"$fastq$R_sufix
    out_p_fw=$OUT_DIR"/02_trimmed_with_sickle/out_paired_"$fastq$F_sufix2
    out_p_rv=$OUT_DIR"/02_trimmed_with_sickle/out_paired_"$fastq$R_sufix2
    out_up_fw=$OUT_DIR"/02_trimmed_with_sickle/out_unpaired_"$fastq$F_sufix2
    out_up_rv=$OUT_DIR"/02_trimmed_with_sickle/out_unpaired_"$fastq$R_sufix2
    out_up=$OUT_DIR"/02_trimmed_with_sickle/out_unpaired_"$fastq$fq_sufix
    echo $fw
    echo $rv

#    sickle pe -t illumina -f $fw -r $rv -o $out_p_fw -p $out_p_rv -s $out_up -q 30 -l 75
    trimmomatic PE $fw $rv $out_p_fw $out_up_fw $out_p_rv $out_up_rv LEADING:3 TRAILING:3 SLIDINGWINDOW:30:30 MINLEN:100
    cat $out_up_fw $out_up_rv > $out_up
    rm $out_up_fw $out_up_rv

    echo "$fastq was trimmed" | tee -a $STD_ERR >> $STD_OUT
done

echo "quality trimming with sickle DONE!"

printf "########## 03_fastqc_trimmed #########\n" | tee -a $STD_ERR >> $STD_OUT
printf "######################################\n" | tee -a $STD_ERR >> $STD_OUT
printf "# fastqc report after trimming\n" | tee -a $STD_ERR >> $STD_OUT
printf "######################################\n" | tee -a $STD_ERR >> $STD_OUT

list=$(\ls $OUT_DIR"/02_trimmed_with_sickle"/* | grep "_"$SLURM_ARRAY_TASK_ID"__")
echo $list
fastqc -o $OUT_DIR/03_fastqc_trimmed -t $SLURM_NTASKS_PER_NODE $list

printf "########## 04_removal_of_16S #########\n" | tee -a $STD_ERR >> $STD_OUT
printf "######################################\n" | tee -a $STD_ERR >> $STD_OUT
printf "# remove 16S reads from trimmed fastq files\n" | tee -a $STD_ERR >> $STD_OUT
printf "######################################\n" | tee -a $STD_ERR >> $STD_OUT
# get input names for bbmap to remove 16_S contamination\n" 2>> $STD_ERR 1>> $STD_OUT
F_sufix="_1.fq"
R_sufix="_2.fq"
# get last character in a string: var="012345" ; echo ${var:0:-1} --> 01234
list=$(for i in `\ls $OUT_DIR"/02_trimmed_with_sickle" | grep out_paired | grep "_"$SLURM_ARRAY_TASK_ID"__"`; do echo ${i:0:-5};done | uniq)
list_up=$(for i in `\ls $OUT_DIR"/02_trimmed_with_sickle" | grep out_unpaired | grep "_"$SLURM_ARRAY_TASK_ID"__"`; do echo ${i:0:-5};done | uniq)

# create output dir and cd into that dir
cd $OUT_DIR/04_removal_of_16S
printf "echo $OUT_DIR"/04_removal_of_16S" # input directory for bbmap
# use bbmap\n" 2>> $STD_ERR 1>> $STD_OUT
for fastq in $list; do
    fw=$OUT_DIR"/02_trimmed_with_sickle/"$fastq$F_sufix
    rv=$OUT_DIR"/02_trimmed_with_sickle/"$fastq$R_sufix
    out_p_fw_mapped=$OUT_DIR"/04_removal_of_16S/mapped_16S_"$fastq$F_sufix
    out_p_rv_mapped=$OUT_DIR"/04_removal_of_16S/mapped_16S_"$fastq$R_sufix
    out_p_fw_unmapped=$OUT_DIR"/04_removal_of_16S/removed_16S_"$fastq$F_sufix
    out_p_rv_unmapped=$OUT_DIR"/04_removal_of_16S/removed_16S_"$fastq$R_sufix

    bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast path=$SILVA_16S_DIR qtrim=rl trimq=10 untrim -Xmx23g in=$fw outm=$out_p_fw_mapped outu=$out_p_fw_unmapped
    bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast path=$SILVA_16S_DIR qtrim=rl trimq=10 untrim -Xmx23g in=$rv outm=$out_p_rv_mapped outu=$out_p_rv_unmapped

done
for fastq in $list_up; do
    up=$OUT_DIR"/02_trimmed_with_sickle/"$fastq$fq_sufix
    out_up_mapped=$OUT_DIR"/04_removal_of_16S/mapped_16S_"$fastq$fq_sufix
    out_up_unmapped=$OUT_DIR"/04_removal_of_16S/removed_16S_"$fastq$fq_sufix
    bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast path=$SILVA_16S_DIR qtrim=rl trimq=10 untrim -Xmx23g in=$up outm=$out_up_mapped outu=$out_up_unmapped
done

printf "########## 05_unaligned_to_hg38 ######\n" | tee -a $STD_ERR >> $STD_OUT
printf "######################################\n" | tee -a $STD_ERR >> $STD_OUT
printf "# remove reads that map to the most \n" | tee -a $STD_ERR >> $STD_OUT
printf "# current human genome release\n" | tee -a $STD_ERR >> $STD_OUT
printf "######################################\n" | tee -a $STD_ERR >> $STD_OUT
# get input names for bbmap to remove 16_S contamination\n" 2>> $STD_ERR 1>> $STD_OUT
F_sufix="_1.fq"
R_sufix="_2.fq"
fq_sufix="__.fq"
# get last character in a string: var="012345" ; echo ${var:0:-1} --> 01234
list=$(for i in `\ls $OUT_DIR"/04_removal_of_16S" | grep removed_16S | grep out_paired | grep "_"$SLURM_ARRAY_TASK_ID"__"`; do echo ${i:0:-5};done | uniq)
list_up=$(for i in `\ls $OUT_DIR"/04_removal_of_16S" | grep removed_16S | grep out_unpaired | grep "_"$SLURM_ARRAY_TASK_ID"__"`; do echo ${i:0:-5};done | uniq)

# create output dir and cd into that dir
cd $OUT_DIR/05_unaligned_to_hg38
echo $OUT_DIR"/05_unaligned_to_hg38" # input directory for bbmap
# use bbmap
for fastq in $list; do
    fw=$OUT_DIR"/04_removal_of_16S/"$fastq$F_sufix
    rv=$OUT_DIR"/04_removal_of_16S/"$fastq$R_sufix
    out_p_fw_unmapped=$OUT_DIR"/05_unaligned_to_hg38/human_clean_"$fastq$F_sufix
    out_p_rv_unmapped=$OUT_DIR"/05_unaligned_to_hg38/human_clean_"$fastq$R_sufix
    out_p_fw_mapped=$OUT_DIR"/05_unaligned_to_hg38/human_al_"$fastq$F_sufix
    out_p_rv_mapped=$OUT_DIR"/05_unaligned_to_hg38/human_al_"$fastq$R_sufix

    bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=$HUMAN_REF_INDEX qtrim=rl trimq=10 untrim -Xmx23g in=$fw outu=$out_p_fw_unmapped outm=$out_p_fw_mapped
    bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=$HUMAN_REF_INDEX qtrim=rl trimq=10 untrim -Xmx23g in=$rv outu=$out_p_rv_unmapped outm=$out_p_rv_mapped

done

for fastq in $list_up; do
    up=$OUT_DIR"/04_removal_of_16S/"$fastq$fq_sufix
    out_up_unmapped=$OUT_DIR"/05_unaligned_to_hg38/human_clean_"$fastq$fq_sufix
    out_up_mapped=$OUT_DIR"/05_unaligned_to_hg38/human_al_"$fastq$fq_sufix
    bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=$HUMAN_REF_INDEX qtrim=rl trimq=10 untrim -Xmx23g in=$up outu=$out_up_unmapped outm=$out_up_mapped
done

printf "###### 06_trimmed_with_sickle_2 ######\n" | tee -a $STD_ERR >> $STD_OUT
printf "######################################\n" | tee -a $STD_ERR >> $STD_OUT
printf "# Make paired end reads the same length\n" | tee -a $STD_ERR >> $STD_OUT
printf "#####################################\n" | tee -a $STD_ERR >> $STD_OUT
printf "# get input names for sickle from removed phiX files\n" | tee -a $STD_ERR >> $STD_OUT
F_sufix="_1.fq"
R_sufix="_2.fq"
fq_sufix="__.fq"
# get last character in a string: var="012345" ; echo ${var:0:-1} --> 01234
list=$(for i in `\ls $OUT_DIR"/05_unaligned_to_hg38" | grep "_"$SLURM_ARRAY_TASK_ID"__"`; do echo ${i:0:-5};done | uniq)
echo 'FASTQ list to be trimmed' | tee -a $STD_ERR >> $STD_OUT
echo $list | tee -a $STD_ERR >> $STD_OUT
echo ""
# create output dir and cd into that dir
cd $OUT_DIR/06_trimmed_with_sickle_2
echo $OUT_DIR"/06_trimmed_with_sickle_2" # input directory for sickle
# use trimmomatic
printf "PAIRED END LOOP TRIMMING STARTED:\n" | tee -a $STD_ERR >> $STD_OUT
for fastq in $list; do
    if [[ "$fastq" == *"_out_paired_"* ]]; then

        fw=$OUT_DIR"/05_unaligned_to_hg38/"$fastq$F_sufix
        rv=$OUT_DIR"/05_unaligned_to_hg38/"$fastq$R_sufix
        out_p_fw=$OUT_DIR"/06_trimmed_with_sickle_2/out_paired_"$fastq$F_sufix
        out_p_rv=$OUT_DIR"/06_trimmed_with_sickle_2/out_paired_"$fastq$R_sufix
        out_up_fw=$OUT_DIR"/06_trimmed_with_sickle_2/out_unpaired_"$fastq$F_sufix
        out_up_rv=$OUT_DIR"/06_trimmed_with_sickle_2/out_unpaired_"$fastq$R_sufix

        #sickle pe -t illumina -f $fw -r $rv -o $out_p_fw -p $out_p_rv -s $out_up_fwrv -q 30 -l 75
        #sickle se -t illumina -f $up -o $out_up -s $out_up -q 30 -l 75
        printf "Trimming paired end read:\n" | tee -a $STD_ERR >> $STD_OUT
        printf "forward: \n$fw\n" | tee -a $STD_ERR >> $STD_OUT
        printf "reverse: \n$rv\n" | tee -a $STD_ERR >> $STD_OUT
        printf "out_p_fw: \n$out_p_fw\n" | tee -a $STD_ERR >> $STD_OUT
        printf "out_p_rv: \n$out_p_rv\n" | tee -a $STD_ERR >> $STD_OUT
        printf "out_up_fw: \n$out_up_fw\n" | tee -a $STD_ERR >> $STD_OUT
        printf "out_up_rv: \n$out_up_rv\n" | tee -a $STD_ERR >> $STD_OUT

        trimmomatic PE $fw $rv $out_p_fw $out_up_fw $out_p_rv $out_up_rv LEADING:3 TRAILING:3 SLIDINGWINDOW:30:30 MINLEN:100
        echo "$fastq was trimmed" | tee -a $STD_ERR >> $STD_OUT
    fi

    if [[ "$fastq" == *"_out_unpaired_"* ]]; then
        up=$OUT_DIR"/05_unaligned_to_hg38/"$fastq$fq_sufix
        out_up=$OUT_DIR"/06_trimmed_with_sickle_2/out_unpaired_"$fastq$fq_sufix
        printf "Trimming unpaired end read:\n" | tee -a $STD_ERR >> $STD_OUT
        printf "up: \n$up\n" | tee -a $STD_ERR >> $STD_OUT
        printf "up_out: \n$up_out\n" | tee -a $STD_ERR >> $STD_OUT
        trimmomatic SE $up $out_up LEADING:3 TRAILING:3 SLIDINGWINDOW:30:30 MINLEN:100
    fi
done
### human clean
up_f_clean=$OUT_DIR"/06_trimmed_with_sickle_2"/$(ls $OUT_DIR"/06_trimmed_with_sickle_2" | grep out_unpaired_human_clean | grep $SLURM_ARRAY_TASK_ID | grep _1.fq)
up_r_clean=$OUT_DIR"/06_trimmed_with_sickle_2"/$(ls $OUT_DIR"/06_trimmed_with_sickle_2" | grep out_unpaired_human_clean | grep $SLURM_ARRAY_TASK_ID | grep _2.fq)
up_clean=$OUT_DIR"/06_trimmed_with_sickle_2"/$(ls $OUT_DIR"/06_trimmed_with_sickle_2" | grep out_unpaired_human_clean | grep $SLURM_ARRAY_TASK_ID | grep __.fq)

cat $up_f_clean $up_r_clean >> $up_clean
rm $up_f_clean $up_r_clean

### human al
up_f_al=$OUT_DIR"/06_trimmed_with_sickle_2"/$(ls $OUT_DIR"/06_trimmed_with_sickle_2" | grep out_unpaired_human_al | grep $SLURM_ARRAY_TASK_ID | grep _1.fq)
up_r_al=$OUT_DIR"/06_trimmed_with_sickle_2"/$(ls $OUT_DIR"/06_trimmed_with_sickle_2" | grep out_unpaired_human_al | grep $SLURM_ARRAY_TASK_ID | grep _2.fq)
up_al=$OUT_DIR"/06_trimmed_with_sickle_2"/$(ls $OUT_DIR"/06_trimmed_with_sickle_2" | grep out_unpaired_human_al | grep $SLURM_ARRAY_TASK_ID | grep __.fq)

cat $up_f_al $up_r_al >> $up_al
rm $up_f_al $up_r_al
#list_up=$(for i in `\ls $OUT_DIR"/05_unaligned_to_hg38" | grep "_"$SLURM_ARRAY_TASK_ID"__" | grep "out_unpaired_"`; do echo ${i:0:-5};done | uniq)
#for fastq in $list_up; do
#    up=$OUT_DIR"/05_unaligned_to_hg38/"$fastq$fq_sufix
#    out_up=$OUT_DIR"/06_trimmed_with_sickle_2/out_unpaired_"$fastq$fq_sufix
#    printf "Trimming unpaired end read:\n" | tee -a $STD_ERR >> $STD_OUT
#    printf $up"\n" | tee -a $STD_ERR >> $STD_OUT
#    echo $up | tee -a $STD_ERR >> $STD_OUT
#
#    trimmomatic SE $up $out_up LEADING:3 TRAILING:3 SLIDINGWINDOW:30:30 MINLEN:100
#    
#    out_up_fw=$OUT_DIR"/06_trimmed_with_sickle_2/out_unpaired_"$fastq$F_sufix
#    out_up_rv=$OUT_DIR"/06_trimmed_with_sickle_2/out_unpaired_"$fastq$R_sufix
#    printf "these files will be combined:\n" | tee -a $STD_ERR >> $STD_OUT
#    echo $out_up_fw
#    echo $out_up_rv
#    cat $out_up_fw $out_up_rv >> $out_up
#    printf "these files will be deleted:\n" | tee -a $STD_ERR >> $STD_OUT
#    echo $out_up_fw
#    echo $out_up_rv
#    rm $out_up_fw $out_up_rv
#done


##############################################################################
### counter end ####
####################
ENDTIME=$(date +%s)
#echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task..."
TOTALTIME=`echo "scale=2;($ENDTIME - $STARTTIME)/60" | bc`  # get time in min by
                                                            # by dividing by 60
echo "It took $TOTALTIME min to complete this task"

