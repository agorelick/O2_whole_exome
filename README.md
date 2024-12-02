# 0. Introduction to the O2 cluster
O2 is a platform for Linux-based high-performance computing at Harvard Medical School. This is the main compute cluster used by the Naxerova Lab at HMS. For reasons of HIPAA compliance and practicality, all “raw” sequencing data should be stored and processed on O2. (Note, “raw” sequencing data refers to sequencing-read-level files like FASTQs, BAMs, SAMs, etc.) Raw sequencing data should never be stored on Dropbox, Google Drive, or other cloud platforms. 

Common operations on raw sequencing data (like mutation calling, copy number profiling, etc.) typically also require large reference data files, which should also be kept on O2. Such reference files may include a version of the human reference genome; and databases of human genetic variation (e.g., common SNPs from gnomAD, dbSNP, 1000 genomes, etc.) As a rule of thumb, all sequencing data that is not processed to the point that you would load into R or open in Excel on a laptop should be worked with on O2. 

To use O2, you will need to register for an account. Instructions are on the O2 wiki: https://harvardmed.atlassian.net/wiki/spaces/O2/overview (See: How to request an O2 account). Requests for new accounts usually take a few days to go through. Note, the wiki for O2 is very well documented and should be your first destination for general information.

## Once you have an O2 account
*	you will need to request that your O2 user be added to the “naxerova” user group. This is necessary for your new user to access the lab’s shared directory on O2 for big data files, which is located at /n/data1/hms/genetics/naxerova/lab/. To request this, submit a ticket on the HMS IT webite (https://harvardmed.service-now.com/stat). _Example ticket description:_ "Hello, I am a postdoc in the Naxerova lab in Department of Genetics at HMS. My account on o2 was recently set up (my username is alg2264). Please add my user to the “naxerova” group, so that I can access our shared data directory at /n/data1/hms/genetics/naxerova/lab. Thank you."
* Create a “scratch” directory (a folder designed for many/large temporarily files that should not be backed up). To create your own scratch directory on O2, please run the following command from a login node (not an interactive job): /n/cluster/bin/scratch_create_directory.sh See: https://harvardmed.atlassian.net/wiki/spaces/O2/pages/2652045313/Scratch+Storage


# 1. QC your FASTQ files 
Before pre-processing your FASTQ files, you can generate a QC report for each file using the commonly used tool FASTQC (www.bioinformatics.babraham.ac.uk/projects/fastqc/). This URL also provides examples of "good" reports and "bad" reports with common problems to aid in your interpretation. 

## Prepare FASTQC slurm script for your data

To run FASTQC on your data, use the `run_fastqc.sh` script included in this directory. Copy this file to a new file for your current data. Make the following edits:
```
## Adjust the run time for your number of fastq files (x2 the number of samples for paired reads). NB: 4 ~1Gb fastq files (running in parallel on 4 cores) takes ~ 5 min to complete.
#SBATCH -t 0-01:00

## replace with your email address
#SBATCH --mail-user=alexander_gorelick@hms.harvard.edu

## change this to the path with FASTQ files to be QC-ed
fq_dir="/n/data1/hms/genetics/naxerova/lab/alex/azenta/30-1098079507_BR3_lpWGS/00_fastq"

## change path for the output files from FASTQC
fastqc_dir="/n/data1/hms/genetics/naxerova/lab/alex/example_patient/fastqc"

## change to your user's scratch directory on O2
tmpdir='/n/scratch/users/a/alg2264/
```

## Run FASTQC as a slurm job
```
# submit the job to slurm
sbatch run_fastqc.sh

# check the status of the job (alg2264 is my O2 username)
squeue -u alg2264
```

## Download/inspect QC results
When FASTQC has finished, copy the resulting directory to your local computer to view the results. From a terminal on your local computer, use `scp` to copy the entire fastqc directory. 
```
# run on your local laptop/desktop
# alg2264 is my username on o2, replace with yours
scp -r alg2264@transfer.rc.hms.harvard.edu:/n/data1/hms/genetics/naxerova/lab/alex/example_patient/fastqc example_fastqc
```
Then you can open the fastqc directory (here, `example_fastqc`) and view a file's results (.html files) in a web browser. For information about how to interpret these results, see: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/


# 2. Run preprocessing pipeline on data

## General workflow
The slurm script `run_preprocessing.sh` provides the steps to preprocess sequencing data using GATK's best practices for variant discovery https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery. I tried to strike a balance between keeping the script simple and also automating as much as possible. The workflow I settled on is to run this script for a single patient at a time, in which all samples for the patient will be processed in parallel. As such, for each new patient, I copy the `run_preprocessing.sh` into a new file and edit the files/information in it for the new specific patient. NB: If later on you have more samples for the same patient, just add their fastq files and sample-names to the file, adjust the `array` argument for the newly added samples, and re-run the script on slurm.

## Re-running run_preprocessing.sh after the job failed
To avoid re-running earlier successful commands in the `run_preprocessing.sh` script, at each step I include logic that checks if the expected output exists before running the command. This way, a sample's job can be resumed mid-way without re-running time-intensive tasks (like bwa). However, this means that if a command fails but generates an output file anyway, the output file must be deleted before that step (and the following steps) can be re-run. If a job fails, make sure you identify where the failure started (in the slurm ".out" output/log file). Then delete the output from that failed job and any subsequent jobs that also failed.

## Edits to make in run_preprocessing.sh for your data
```#!/bin/bash
#SBATCH -c 20           # number of cores requested for each sample. Can be changed, but it should align with all core numbers specified in the commands below
#SBATCH -t 0-4:00       # ~24-36hrs for each Exome sample. ~12hrs for lpWGS samples (~1x)
#SBATCH -p short        # > 12 hrs requires "medium" queue, <= 12 hrs can be on "short" queue
#SBATCH --mem 16000     # memory requested for each sample in Mb
#SBATCH --array=0-2     # this should be a 0-indexed array with a length equal to the number of samples. Here: "0,1,2" for 3 samples
#SBATCH -o run_preprocessing_%A_%a.out
#SBATCH --mail-user=alexander_gorelick@hms.harvard.edu  # Replace with your email address. Slurm notifications will be emailed here
#SBATCH --mail-type=ALL

# paths for output files.
cutadapt_path='/n/data1/hms/genetics/naxerova/lab/alex/example_patient/cutadapt'
preprocessing_path='/n/data1/hms/genetics/naxerova/lab/alex/example_patient/preprocessing'
bam_path='/n/data1/hms/genetics/naxerova/lab/alex/example_patient/bams'

# path to the human reference genome. Here, version GRCh38 (aka hg38).
# To use with bwa mem (for alignment of fastq files), this must first be indexed using BWAIndex/version0.6.0 (needs to be greater than 0.5.* for bwa>=0.7)
reference='/n/data1/hms/genetics/naxerova/lab/alex/reference_data/assemblies/Homo_sapiens_NCBI_GRCh38/NCBI/GRCh38/Sequence/BWAIndex/genome.fa'

# path to a database of common SNPs for the appropriate reference genome version (required for
polymorphic_sites='/n/data1/hms/genetics/naxerova/lab/alex/reference_data/dbSNP/dbSNP_GRCh38/00-common_all_renamedchrs.vcf.gz'

# path to a directory where many large temporary files can safetly be written. On O2, this is best done in your "scratch" directory (see O2 wiki for instructions to create this directory).
tmpdir='/n/scratch/users/a/alg2264/'

# adapter sequences based on library prep.
# the sequences here are for "Illumina universal adapter sequences" used for both Twist Human Core Exome kits and Illumina TruSeq kits (the WES and lpWGS kits used by Azenta respectively)
# If you aren't sure what adapter sequence was used, check the output of fastqc.
adapter_seq_R1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
adapter_seq_R2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'

# Original FASTQ file suffices+extensions
suffix_R1='_R1_001.fastq.gz'
suffix_R2='_R2_001.fastq.gz'

# update this array with the "prefix" for each sample's FASTQ data files.
# Example: if you have paired-end FASTQ files "path/s1_R1_001.fastq.gz" and "path/s1_R2_001.fastq.gz" for a sample "s1", the prefix should be "path/s1"
declare -a fq_prefices=(
    "/n/data1/hms/genetics/naxerova/lab/alex/azenta/30-1098079507_BR3_lpWGS/00_fastq/1A"
    "/n/data1/hms/genetics/naxerova/lab/alex/azenta/30-1098079507_BR3_lpWGS/00_fastq/1B"
    "/n/data1/hms/genetics/naxerova/lab/alex/azenta/30-1098079507_BR3_lpWGS/00_fastq/3B")

# update this array with the name for each sample in the order of the FASTQ prefices.
declare -a samples=("Liv1a-A" "Liv1b-A" "N1")
```

## Run preprocessing pipeline as a slurm job
```
# submit the job to slurm
sbatch run_preprocessing.sh

# check the status of the job (alg2264 is my O2 username)
squeue -u alg2264
```

