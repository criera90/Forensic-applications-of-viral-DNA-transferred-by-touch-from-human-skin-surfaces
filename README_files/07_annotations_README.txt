All the final taxonomy files are within the folder /PATH_TO_YOUR_WORKING_DIR/07_annotations/07_summaries. Below is a guide on how to replicate the process:

1. Copy the "unmapped_TO_clean_assemblies_100_negs_gt100bp_CONTIGS.fa" file from `PATH_TO_YOUR_WORKING_DIR/02_clean_assemblies/output/03_removal_of_neg_control_contigs` to every input folder within all these directories:

/PATH_TO_YOUR_WORKING_DIR/07_annotations/07_annotation_checkV        
/PATH_TO_YOUR_WORKING_DIR/07_annotations/07_annotations_blast        
/PATH_TO_YOUR_WORKING_DIR/07_annotations/07_annotations_kraken       
/PATH_TO_YOUR_WORKING_DIR/07_annotations/07_annotations_mapping_to_db
/PATH_TO_YOUR_WORKING_DIR/07_annotations/07_annotations_metaphlan    
/PATH_TO_YOUR_WORKING_DIR/07_annotations/07_annotations_vibrant      
/PATH_TO_YOUR_WORKING_DIR/07_annotations/07_summaries

2. Download the NCBI Virus Assembly database in fasta format. Copy the file to `PATH_TO_YOUR_WORKING_DIR/07_annotations/07_annotations_mapping_to_db/input`

3. Run the scripts inside the scrip folder for each of the annotation tools:

/PATH_TO_YOUR_WORKING_DIR/07_annotations/07_annotation_checkV/scripts        
/PATH_TO_YOUR_WORKING_DIR/07_annotations/07_annotations_blast/scripts        
/PATH_TO_YOUR_WORKING_DIR/07_annotations/07_annotations_kraken/scripts       
/PATH_TO_YOUR_WORKING_DIR/07_annotations/07_annotations_mapping_to_db/scripts
/PATH_TO_YOUR_WORKING_DIR/07_annotations/07_annotations_metaphlan/scripts    
/PATH_TO_YOUR_WORKING_DIR/07_annotations/07_annotations_vibrant/scripts      
                 

4. Change directory to `/PATH_TO_YOUR_WORKING_DIR/07_annotations/07_summaries` and follow the following instructions to unify annotated contigs:

Download the taxonomy files from the NCBI Taxonomy FTP site: https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/
Copy the folder to any location in your working directory. You'll need the path to this directory for the scripts using the Python taxidtools (https://github.com/CVUA-RRW/taxidtools)

Kraken
Taken directly from the kraken_output.txt, which includes taxids. Used taxidtools to get full taxonomy. (https://github.com/CVUA-RRW/taxidtools)

cat ../input/kraken_output.txt | grep -v Homo | grep vir | awk -F "\t" '{print $2"\t"$3}' | sed 's/(taxid /\t/g' | awk -F ")" '{print $1}' > kraken_contigs__contig_species_taxid__.txt

Vibrant
Taken from vibrant output. Used the ones that were not already classified by kraken. Used blast to get taxid and taxidtools to get full taxonomy.

cat ../input/VIBRANT_machine_unmapped_TO_clean_assemblies_100_negs_gt100bp_CONTIGS.tsv | grep virus | awk -F "\t" '{print $1}' > vibrant_contigs.txt

cat kraken_contigs.txt | sort > kraken_contigs_sorted.txt
cat vibrant_contigs.txt | sort > vibrant_contigs_sorted.txt
comm -13 kraken_contigs_sorted.txt vibrant_contigs_sorted.txt | sort > only_vibrant_contigs_sorted.txt

# contigs classified by vibrant but not by QMAP. Because vibrant doesn’t provide tax, and qmap assigns acc. Assigning taxonomy just by blast (method for vibrant), is less accurate than mapping directly to the database. 

comm -23 only_vibrant_contigs_sorted.txt qmap_contigs_sorted.txt | sort > only_vibrant_contigs_NOT_IN_QMAP_sorted.txt

Checkv
Taken from quality report of checV. Used the ones with at least 1 viral gene and that were not classified by kraken or vibrant. Used blast to get taxid and taxidtools to get full taxonomy.

cat ../input/quality_summary.tsv | awk -F "\t" '{ if ($6 >= 1) print $1}' > checkV_contigs.txt

comm -13 kraken_vibrant_QMAP_final_sorted.txt checkV_contigs_sorted.txt | sort > checkV_NOT_IN_KRAKEN_VIBRANT_QMAP_sorted.txt

QMAP
Used those with QMAP higher than or equal to 20 and that were not classified by kraken, vibrant or checkv. 

cat ../input/very_sensitive_local_unmapped_TO_clean_assemblies_100_negs_gt100bp_CONTIGS.fa.sam | awk -F "\t" '{if ($5 >=15) print $1"\t"$3"\t"$5}' > qmap_contigs__contig_acc_qmap__.txt

comm -13 kraken_contigs_sorted.txt qmap_contigs_sorted.txt | sort > qmap_contigs_NOT_IN_KRAKEN_sorted.txt

cat kraken_vibrant_QMAP_final_sorted.txt checkV_NOT_IN_KRAKEN_VIBRANT_QMAP_sorted.txt | sort > classified_sorted.txt

The following logic was used to select taxonomy assignemnts from each tool. 
	
kraken = 1
kraken – vibrant = 2
qmap – vibrant = 3
1 + 2 + 3 = 4
checkV – 4 = 5

1+2+3+5 = classified

Blast_passed - Classified = 6
comm -13 classified_sorted.txt blast_passed_sorted.txt | sort > blast_passed2_sorted.txt

# esearch -db nucleotide -query "NC_012783.2" | esummary | xtract -pattern DocumentSummary -element TaxId
# epost -db nucleotide -id NC_012783.2, NC_022098.1 | efetch -format uid | esummary -db nucleotide | xtract -pattern DocumentSummary -element TaxId,AccessionVersion
# get taxids with this command (must have entrez-direct installed)
# list of e-direct database names: https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
# for i in `echo $b`;do tax=$(esearch -db nucleotide -query $i | esummary | xtract -pattern DocumentSummary -element TaxId); printf "$i\t$tax\n" >> seqid_to_taxid.txt;done


python code to get contigs from qmap that are not present in kraken
First need to get sam file with the first 10 fields (from the 07_annotations_mapping_to_db output. This is because the 11th field has symbols that mess up with Pandas' import functions. 

cat very_sensitive_local_unmapped_TO_clean_assemblies_100_negs_gt100bp_CONTIGS.fa.sam | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > very_sensitive_local_unmapped_TO_clean_assemblies_100_negs_gt100bp_CONTIGS_short.fa.sam

import pandas as pd
qmap_dir="very_sensitive_local_unmapped_TO_clean_assemblies_100_negs_gt100bp_CONTIGS.fa.sam"
qmap_contigs_dir="qmap_contigs_NOT_IN_KRAKEN_sorted.txt"
col_names=["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ"]
qmap=pd.read_csv(qmap_dir,sep="\t", header=None, names=col_names)
qmap_contigs=pd_csv(qmap_contigs_dir, sep="\t", header=None, names=["contigs"])
qmap[qmap["QNAME"].isin(qmap_contigs["contigs"])][["QNAME","RNAME"]].to_csv("qmap__QNAME_RNAME__.txt",sep="\t",index=False)

On the command line
b=$(cat qmap__QNAME_RNAME__2.txt | awk -F "\t" '{print $2}' | sed ':a;N;$!ba;s/\n/ /g')
module load entrez-direct/16.2
for i in `echo $b`;do tax=$(esearch -db nucleotide -query $i | esummary | xtract -pattern DocumentSummary -element TaxId); printf "$i\t$tax\n" >> qmap__RNAME__TAXID.txt;done

paste qmap__QNAME_RNAME__2.txt qmap__RNAME__TAXID.txt | awk -F "\t" '{print $1"\t"$2"\t"$4}' > qmap__QNAME_RNAME__TAXID.txt

Format correctly
cut --complement -d',' -f2,3 ../input/tax_checkV.txt > tax_checkV_ready.txt
cut --complement -d',' -f2 ../input/tax_kraken.txt > tax_kraken_ready.txt
cut --complement -d',' -f2,3 ../input/tax_qmap.txt > tax_qmap_ready.txt
cut --complement -d',' -f2,3 ../input/tax_remaining.txt > tax_remaining_ready.txt
cut --complement -d',' -f2,3 ../input/tax_vibrant2.txt > tax_vibrant_ready.txt

Remaining. 
Unclassified contigs were assigned taxonomy as follows. Annotation of remaining contigs that were not annotated with kraken, vibrant, checkV or mapping to NCBI viral database with bowtie2: The list of all contigs was pulled from the ‘kraken_output.txt’ file. Only those classified as non-human by kraken and not classified as virus by other annotation tool (kraken, vibrant, checkV and mapping to NCBI viral database) were selected to annotate with BLAST (evalue < 1e-17, piden > 95% and qcov > 95%). 

cat ../input/kraken_output.txt | grep U | awk -F "\t" '{print $2}' > unclasified_by_kraken.txt
cat TO_blast.txt | awk -F "\t" '{if ($5 <= 1e-17 && $6 > 95 && $7 > 95) print $1}' | sort | uniq | wc -l


PATH_TO_YOUR_WORKING_DIR
