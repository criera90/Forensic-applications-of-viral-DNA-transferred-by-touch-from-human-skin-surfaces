1. Download the datasets from NCBI from the BIOPROJECT PRJNA1220097

2. Copy all the file in PATH_TO_YOUR_WORKING_DIR/00_QC/input/

3. Modify the file names so that each forward and reverse file has an index (from 1 to total number of samples) preceding the filename and two underscores. It should look like this:

.
.
.
19__TO_10_1_Ceramic_1.fq.gz
19__TO_10_1_Ceramic_2.fq.gz
20__TO_10_1_Glass_1.fq.gz
20__TO_10_1_Glass_2.fq.gz
21__TO_10_1_Hands_1.fq.gz
21__TO_10_1_Hands_2.fq.gz
22__TO_10_1_Knife_1.fq.gz
22__TO_10_1_Knife_2.fq.gz
23__TO_10_1_Neg_1.fq.gz
23__TO_10_1_Neg_2.fq.gz
24__TO_10_1_P10D1_1.fq.gz
24__TO_10_1_P10D1_2.fq.gz
.
.
.

The files name format has the following (using 19__TO_10_1_Ceramic_1.fq.gz as example):
Index		"19"
Separator 	"__" (two underscores)
Sample type	"TO" (for touch objects)
Separator 	"_" (one underscore)
Person_ID	"10"
Separator 	"_" (one underscore)
Day of collection	"1"
Separator 	"_" (one underscore)
Object type	"Ceramic"
Separator 	"_" (one underscore)
Fw/Rv sequence	"1" or "2" for forward of reverse Illumina pair
Separator 	"."
file extension	"fq.gz" (compressed fastq file)


4. Modify the script 'mg_quality_check_TO_redo.sh' to match your workspace input paths and filenames. Also might need to modify memory and the number of cores in the SLURM script

5. The final output should be free of host contamination and low-quality reads. 

