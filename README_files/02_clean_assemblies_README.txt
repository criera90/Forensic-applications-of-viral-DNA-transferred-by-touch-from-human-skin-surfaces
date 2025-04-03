1. Copy the final samples co-assembly file named "final.contigs.fa" from `PATH_TO_YOUR_WORKING_DIR/01_assemblies/output/TO_pooled_samples/megahit_out` to `PATH_TO_YOUR_WORKING_DIR/02_clean_assemblies/input` and rename it "TO_pooled_assembly_samples.fa"

2. 1. Copy the final negative samples co-assembly file named "final.contigs.fa" from `PATH_TO_YOUR_WORKING_DIR/01_assemblies/output/TO_pooled_negs/megahit_out` to `PATH_TO_YOUR_WORKING_DIR/02_clean_assemblies/input` and rename it "TO_pooled_assembly_negs.fa"

3. Modify the script 'mg_TO_clean_assemblies.sh' to match your workspace input paths and filenames. Also might need to modify memory and the number of cores in the SLURM script

5. The final output should be a fasta file with contigs not present in negative controls called "unmapped_TO_clean_assemblies_100_negs_gt100bp_CONTIGS.fa"



