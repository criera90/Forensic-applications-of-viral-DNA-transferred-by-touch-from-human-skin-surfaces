1. Copy bam files generated in the previous step from `PATH_TO_YOUR_WORKING_DIR/03_mapping_reads_back_to_contigs/output/01_mapped_contigs` to `PATH_TO_YOUR_WORKING_DIR/04_getting_count_table/input`

2. Run the script "mg_get_count_table_TO.sh" modifying pathnames as necessary.

3. Change directory to `PATH_TO_YOUR_WORKING_DIR/04_getting_count_table/output` and run the script (make sure you are using Python 2.7):

module load python/2.7

python ../../scripts/get_count_table.py 01_idxstats/* > TO_counts.txt
sed 's/[[:digit:]]\+\__//g' TO_counts.txt | sed 's/.sorted.bam.idxstats.txt//g' | sed 's/01_idxstats\///g' | sed 's/Ceram/Ceramic/g' > TO_counts_edited.txt

4. File generated "TO_counts_edited.txt" is the count table for the subsequent analyses using phyloseq R package. 