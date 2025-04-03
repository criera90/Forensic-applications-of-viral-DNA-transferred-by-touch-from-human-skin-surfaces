1. Copy the sample files (Ceramic, Glass, Plastic, Knife, Phone) having "_human_clean_*" in the name from `PATH_TO_YOUR_WORKING_DIR/00_QC/output/06_trimmed_with_sickle_2` to `PATH_TO_YOUR_WORKING_DIR/01_assemblies/input/TO_samples`

2. Copy the negative sample files having "_human_clean_*" in the name from `PATH_TO_YOUR_WORKING_DIR/00_QC/output/06_trimmed_with_sickle_2` to `PATH_TO_YOUR_WORKING_DIR/01_assemblies/input/TO_negs`

3. Modify the script 'mg_assemble_with_megahit_TO_pooled_negs.sh' to match your workspace input paths and filenames. Also might need to modify memory and the number of cores in the SLURM script

4. Modify the script 'mg_assemble_with_megahit_TO_pooled_samples.sh' to match your workspace input paths and filenames. Also might need to modify memory and the number of cores in the SLURM script

5. The final output should be a co-assembly of all samples and a co-assembly of all negative samples

