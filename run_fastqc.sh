#!/bin/bash
#SBATCH -c 4                # number of cores requested for FASTQC. This number of fastq files will be run in parallel
#SBATCH -t 0-01:00          # total run time for all fastq files combined. "0-01:00" means 1 hr. 4 ~1Gb fastq files (running in parallel on 4 cores) takes ~ 5 min.
#SBATCH -p short            # > 12 hrs requires "medium" queue, <= 12 hrs can be on "short" queue
#SBATCH --mem 2000          # memory requested in Mb
#SBATCH -o run_fastqc.out   # slurm log file
#SBATCH --mail-user=alexander_gorelick@hms.harvard.edu  # Email to which notifications will be sent
#SBATCH --mail-type=ALL

# load pre-installed FASTQC module on O2 
module load fastqc/0.11.9

# path to FASTQ files to be QC-ed
fq_dir="/n/data1/hms/genetics/naxerova/lab/alex/azenta/30-1098079507_BR3_lpWGS/00_fastq"

# path to output files from FASTQC
fastqc_dir="/n/data1/hms/genetics/naxerova/lab/alex/example_patient/fastqc"

# path to a directory where many large temporary files can safetly be written. On O2, this is best done in your "scratch" directory (see O2 wiki for instructions to create this directory).
tmpdir='/n/scratch/users/a/alg2264/'

# create directory for FASTQC output if it doesn't exist
mkdir -p $fastqc_dir

# run FASTQC. Edit the file extension to match your fastq files (here, ".fastq.gz")
fastqc ${fq_dir}/*.fastq.gz --threads 4 --dir $tmpdir --outdir $fastqc_dir


