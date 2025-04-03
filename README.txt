This repository contains instructions for replicating the data analysis of "Forensic applications of viral DNA transferred by touch from human skin surfaces."

The first step is to git clone the repository in your working directory. 

The steps for processing the raw sequences are described in the README files in the README_files folder. 

The repository also structures each step of the analysis in folders enumerated from 00 to 07. Before running the scripts, you must `mkdir input output scripts SLURM_reports` folders within each folder. Not doing this would create problems given the way the scripts are coded.  

The "scripts" folders within each directory contain the script that processes input files and outputs intermediary files for the next steps. Each step has a README file in the README_files folder explaining how to download and what files need to be copied from the output file of the previous steps into the input folder of the next step. 

The final annotation files are in the 07_summaries folder within the 07_annotations folder, except the blast outputs due to the file size. The Blast output should be obtained following the instructions in the Blast annotation step. 
The quality_summary.tsv from checkV and the kraken_output.txt files in `/Users/crieraruiz/Documents/submit_scritps_to_github/07_annotations/07_summaries/00_get_viral_contigs_classified_by_tools` where not included due to file size and will have to be copied from the output of each tool after you run it. 

The consolidated files required to replicate the figures are in the "Paper figures" folder. R is required to create each figure. 

If there's any question, please refer to the original repository (), or contact the first author:

Carlos Riera-Ruiz
criera-ruiz@huskers.unl.edu
carlos.riera.ruiz@gmail.com
criera@espol.edu.ec