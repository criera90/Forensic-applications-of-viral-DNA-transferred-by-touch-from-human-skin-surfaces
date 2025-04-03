1. Download the most recent version of the human genome: https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=9606

2. Download the most recent SSU dataset from te SILVA database: https://www.arb-silva.de/download/arb-files

or

# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz .
# wget https://ftp.arb-silva.de/release_138.1/Exports/SILVA_138.1_SSURef_tax_silva.fasta.gz .

3. Copy both file in ~/00_indices/input/human and ~/00_indices/input/silva, respectively

4. Modify the script \'91mg_make_16S_and_human_index_with_bbmap.sh\'92 to match your workspace\'92s input paths and filenames

