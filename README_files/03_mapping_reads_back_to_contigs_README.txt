1. Copy the clean sample files from `PATH_TO_YOUR_WORKING_DIR/00_QC/output/06_trimmed_with_sickle_2` to `03_mapping_reads_back_to_contigs/input/TO_samples`

2. remove the "out_unpaired_human_clean_removed_16S_out_unpaired_phiX_removed_" or "out_paired_human_clean_removed_16S_out_unpaired_phiX_removed_" string from the file names. The files now should look like this:

15__TO_10_1_P10D1_2.fq  
15__TO_10_1_P10D1__.fq  
16__TO_10_1_Phone_1.fq  
16__TO_10_1_Phone_2.fq  
16__TO_10_1_Phone__.fq  
17__TO_10_1_Plastic~.fq 
17__TO_10_1_Plastic~.fq 
17__TO_10_1_Plastic__.fq
.
.
.

The ones ending with 1.fq are the forward
The ones ending with 2.fq are the reverse
The ones ending with __.fq are the unpaired sequences after quality trimming. 

3. Copy the "unmapped_TO_clean_assemblies_100_negs_gt100bp_CONTIGS.fa" file from "PATH_TO_YOUR_WORKING_DIR/02_clean_assemblies/output/03_removal_of_neg_control_contigs"

4. First, run the "mg_index_mapping_contigs_TO.sh" script and modify path directories as necessary to match your work directory. 

5. Then run the "mg_mapping_contigs_TO.sh" modifying pathnames as necessary.


5. The final output should be a fasta file with contigs not present in negative controls called "unmapped_TO_clean_assemblies_100_negs_gt100bp_CONTIGS.fa"



