#!/bin/sh
#SBATCH --job-name=bay_scallop_filter_test # Job name
#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=d.ence@mail.ufl.edu  # Where to send mail	
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=1                   # Run a single task	
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=10gb                  # Total memory limit
#SBATCH --time=96:00:00              # Time limit hrs:min:sec
#SBATCH --out=bay_scallop_filter_%j.out     # Standard output and error log
#SBATCH --err=bay_scallop_filter_%j.err     # Standard output and error log
#SBATCH --qos=kirst-b
#SBATCH --partition=hpg1-compute

#export OMP_NUM_THREADS=4
 
# Load required modules; for example, if your program was
# compiled with Intel compiler, use the following 
module load python
module load fastqc
 
# Run your program with correct path and command line options

FILTER_SCRIPT="/home/d.ence/applications/10X_fastq_filter/FastqFilter_10X.py"
BAD_SEQ_FILE="/home/d.ence/applications/10X_fastq_filter/test_data/R2badseqs.txt"

FASTQ_1="5452-MK-0001_S4_L003.filtered.compl_test.1.fastq.gz"
FASTQ_2="5452-MK-0001_S4_L003.filtered.compl_test.2.fastq.gz"
fastqc -f fastq $FASTQ_1 $FASTQ_2

#python $FILTER_SCRIPT --gzip -1 $FASTQ_1 -2 $FASTQ_2 -f $BAD_SEQ_FILE -p 2 -o 5452-MK-0001_S4_L003.filtered.aa.compl_test

#FILTERED_FASTQ_1="5452-MK-0001_S4_L003.filtered.aa.compl_test.1.fastq.gz"
#FILTERED_FASTQ_2="5452-MK-0001_S4_L003.filtered.aa.compl_test.2.fastq.gz"

#fastqc -f fastq $FILTERED_FASTQ_1 $FILTERED_FASTQ_2 

