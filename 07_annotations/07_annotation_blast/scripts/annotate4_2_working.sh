#!/bin/sh

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
#                   -word_size 28 \
#                   -evalue 0.5 \
#                   -reward 1 \
#                   -penalty -2 \
#                   -gapopen 0 \
#                   -gapextend 2.5 \     
######################################################################################
######################################################################################
# counter start
STARTTIME=$(date +%s)

# https://www.baeldung.com/linux/use-command-line-arguments-in-bash-script
# try this one to fix -h output: 
# https://tecadmin.net/pass-command-line-arguments-in-shell-script/
while getopts i:f:o:h flag
do
    case "${flag}" in
        i) file=${OPTARG};;
        f) blast_file=${OPTARG};;
        o) out_path=${OPTARG};;
        h) echo "
            Usage:
            sbatch annotate.sh -i table -f fasta_file_with_actual_seqs

            -i      address of a fasta file containing the sequences of the queries
                    if the fasta file is in the same directory where the program is
                    being excecuted, no need to add directory

            -f      fasta file with sequences of acc in input

            -h      help

            Notes:
            sufix   defined as a keyword that you want to use to filter the hits 
                    on your description table 
                    (at the moment the field used by default is Scientific name)
                    "
          exit;;
    esac
done

#cd .
#echo $PWD
#module load blast/2.10
#module load biodata
#export BLASTDB=/work/HCC/BCRF/BLAST     # this line is important to make sure
                                        # BLAST will know where to locate the 
                                        # taxonomy database. This is a problem 
                                        # specially when we ask for sscinames. Blast
                                        # will look only into the nt if $BLAST/nt is 
                                        # specified, but not into the taxbd, which is
                                        # also in the $BLAST directory

out_name=$(echo $blast_file | awk -F "/" '{print $NF}' | cut -d. -f1)

# BLAST fasta_file list against nt database
#blastn -query $fasta_file \
#    -out $out_path/$out_name"_blast.txt" \
#    -db $BLAST/nt \
#    -outfmt "6 qaccver saccver sallseqid sscinames evalue pident qcovs staxids sskingdoms" \
#    -num_threads 16 \
#    -max_hsps 200 \
#    -word_size 28 \
#    -evalue 10 \
#    -reward 1 \
#    -penalty -2 \
#    -gapopen 0 \
#    -gapextend 2

# try https://cran.r-project.org/web/packages/taxonomizr/readme/README.html

#while [[ ! -f $out_name"_blast.txt" ]]; do
#    sleep 1
#done
#
#
# parallelize the job:
# https://unix.stackexchange.com/questions/103920/parallelize-a-bash-for-loop
# how to use sed in the scinames last pipe:
# https://linuxhint.com/newline_replace_sed/
#for table in `\ls $file`; do
while read -r i; do
    STARTTIME_=$(date +%s)
    subject=$(echo "$i" | awk -F " " '{print $1}')
    subject=$(echo "$i")
    echo "blast file is $blast_file"
    printf "$subject"
    scinames=$(cat $blast_file | awk -v var="$subject" -F "\t" '{if ($1 == var && $5 < 0.0000000001 && $6 > 95 && $7 > 95) print $4 "\t" $8}' | head -1 | uniq | sed ':a;N;$!ba;s/\n/,/g')
#    scinames=$(cat $blast_file | awk -v var="$subject" -F "\t" '{if ($1 == var) print $4}' | head -3 | sed ':a;N;$!ba;s/\n/,/g')
    #scinames=$(cat $blast_file | grep $subject)
    echo "These are the scinames $scinames"
    printf "$i\t$scinames\n" >> annotated_best_hit_sscinames_staxids.txt
    ENDTIME_=$(date +%s)
done < $file
#done
# counter end
ENDTIME=$(date +%s)
#echo "It takes $(($ENDTIME - $STARTTIME)) seconds to complete this task..."
TOTALTIME=`echo "scale=2;($ENDTIME - $STARTTIME)/60" | bc`  # get time in min by
                                                            # by dividing by 60
echo "It took $TOTALTIME min to complete this task"

# notes:
# to use variables inside sed commands, you need to use double quotes
# (https://askubuntu.com/questions/76808/how-do-i-use-variables-in-a-sed-command)

