#!/bin/sh

#SBATCH --time=168:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=10GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --job-name=faster_TO_blast_taxonomy
#SBATCH --error=/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_blast/SLURM_reports/job.%J.%x.err
#SBATCH --output=/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_blast/SLURM_reports/job.%J.%x.out


######################################################################################
######################################################################################
# Program name      annotate3.sh
# Version           v1.3
# Date modified     211201
# Function          outputs the description table from a BLASTn search
#                   uses as input a file with fasta formatted sequences
#                   the output is formatted in outfmt 6 with the following fields:
#                   qaccver         Query accession version
#                   saccver         Subject accessio version
#                   sallseqid       All subject Seq-id(s), separated by a ';'
#                   salltitles      All Subject Title(s), separated by a '<>'
#                   evalue          Expect value
#                   pident          Percentage of identical matches
#                   qcov            query coverage
#                   ssciname        subject scientific name
#
# Notes:            this version uses a normal tblastn parameters (highly similar
#                   sequences)
#                   -word_size 11 \
#                   -evalue 0.5 \
#                   -reward 2 \
#                   -penalty -3 \
#                   -gapopen 5 \
#                   -gapextend 2 \     
#                   
#                   To get the taxids lists to limit search to viruses:
#                   Limiting a BLAST search with a high-level taxonomic node
#                   
#                   $ get_species_taxids.sh -n Enterobacterales
#                   Taxid: 91347
#                    rank: order
#                    division: enterobacteria
#                    scientific name: Enterobacterales
#                    common name:
#                   1 matches found
#                   $ get_species_taxids.sh -t 91347 > 91347.txids
#                   $ blastn –db nt –query QUERY –taxidlist 91347.txids –outfmt 7 –out OUTPUT.tab
#######################################################################################
######################################################################################

# module load
module load blast/2.14

# defining directories
INPUT_DIR="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_blast/input"
OUT_DIR="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2_TO_2/07_annotations/07_annotations_blast/output/faster"
fasta_file="unmapped_TO_clean_assemblies_100_negs_gt100bp_CONTIGS.fa"
txids="viruses_txids.txt"
out_name=$(echo $fasta_file | awk -F "/" '{print $NF}' | cut -d. -f1)
db_dir="/lustre/work/adamowiczlab/crieraruiz/PROJECTS/NIJ2/07_annotations_gt1000/07_gt1000_blast/input/NCBI"

STARTTIME=$(date +%s)
# copy db to scratch
cp $db_dir/nt_viruses* /scratch
# counter end
ENDTIME=$(date +%s)
#echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task..."
TOTALTIME=`echo "scale=2;($ENDTIME - $STARTTIME)/60" | bc`  # get time in min by

echo "It took $TOTALTIME min to complete the db copy to scratch"

# counter start
STARTTIME=$(date +%s)

#    -taxidlist $INPUT_DIR/$txids \
#    -mt_mode 1 \
# BLAST fasta_file list against nt database
blastn -query $INPUT_DIR/$fasta_file \
    -out $out_name"_blast.txt" \
    -db /scratch/nt_viruses \
    -outfmt "6 qaccver saccver sallseqid sscinames evalue pident qcovs staxids sskingdoms" \
    -num_threads $SLURM_NTASKS_PER_NODE \
    -max_hsps 1 \
    -word_size 11 \
    -evalue 1 \
    -reward 2 \
    -penalty -3 \
    -gapopen 5 \
    -gapextend 2

cp $out_name"_blast.txt" $OUT_DIR/

#while [[ ! -f $out_name"_blast.txt" ]]; do
#    sleep 1
#done
#
#
## parallelize the job:
## https://unix.stackexchange.com/questions/103920/parallelize-a-bash-for-loop
## how to use sed in the scinames last pipe:
## https://linuxhint.com/newline_replace_sed/
#for table in `\ls $file`; do
#    while read -r i; do
#        STARTTIME_=$(date +%s)
#        subject=$(echo "$i" | awk -F " " '{print $1}')
#        echo "blast file is $out_name_blast.txt"
#        echo "subject is $subject"
#        scinames=$(cat $out_name"_blast.txt" | grep $subject | 
#            awk -F "\t" '{if ($5 == 0 && $6 == 100 && $7 == 100); {print $4}}' | sed ':a;N;$!ba;s/\n/,/g')
#        echo "This are the scinames $scinames"
#        printf "$i\t$scinames\n" >> $table"_annotated.txt"
#        ENDTIME_=$(date +%s)
#    done < $file/$table
#done

# counter end
ENDTIME=$(date +%s)
#echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task..."
TOTALTIME=`echo "scale=2;($ENDTIME - $STARTTIME)/60" | bc`  # get time in min by
                                                            # by dividing by 60
echo "It took $TOTALTIME min to complete this BLAST"

# notes:
# to use variables inside sed commands, you need to use double quotes
# (https://askubuntu.com/questions/76808/how-do-i-use-variables-in-a-sed-command)

