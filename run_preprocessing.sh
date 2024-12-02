#!/bin/bash
#SBATCH -c 20           # number of cores requested for each sample. Can be changed, but it should align with all core numbers specified in the commands below
#SBATCH -t 0-4:00       # ~24-36hrs for each Exome sample. ~12hrs for lpWGS samples (~1x)
#SBATCH -p short        # > 12 hrs requires "medium" queue, <= 12 hrs can be on "short" queue
#SBATCH --mem 16000     # memory requested for each sample in Mb
#SBATCH --array=0-2     # this should be a 0-indexed array with a length equal to the number of samples. Here: 0,1,2 for 3 samples
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

# update this array with the "prefix" for each sample's FASTQ data files. 
# Example: if you have paired-end FASTQ files "path/s1_R1_001.fastq.gz" and "path/s1_R2_001.fastq.gz" for a sample "s1", the prefix should be "path/s1"
declare -a fq_prefices=(
    "/n/data1/hms/genetics/naxerova/lab/alex/azenta/30-1098079507_BR3_lpWGS/00_fastq/1A"
    "/n/data1/hms/genetics/naxerova/lab/alex/azenta/30-1098079507_BR3_lpWGS/00_fastq/1B"
    "/n/data1/hms/genetics/naxerova/lab/alex/azenta/30-1098079507_BR3_lpWGS/00_fastq/3B")

# update this array with the name for each sample in the order of the FASTQ prefices.
declare -a samples=("Liv1a-A" "Liv1b-A" "N1")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define paths for this sample's input/output files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# the path and sample-name for the current sample (this lets us submit jobs for each sample in parallel by using the "--array" argument in the slurm header)
# let's also define sam-file read groups/tags/meta-data here as some are sample-specific 
i=${SLURM_ARRAY_TASK_ID}
fq_prefix=${fq_prefices[$i]}
sample=${samples[$i]}   # sample name, will be used as the "SM" value in the resulting bam file
readgroup=$i            # read-tag for sample
PL='Illumina'           # read-tag for platform/technology used

# specify intermediate/output files for this sample
fq_R1="${fq_prefix}_R1_001.fastq.gz"
fq_R2="${fq_prefix}_R2_001.fastq.gz"
cutadapt_R1="${cutadapt_path}/trimmed.${sample}_R1_001.fastq"
cutadapt_R2="${cutadapt_path}/trimmed.${sample}_R2_001.fastq"
sam_file="${preprocessing_path}/${sample}_raw.sam"
marked_dup_bam_file="${preprocessing_path}/${sample}_marked_dup.bam"
marked_dup_metrics_file="${preprocessing_path}/${sample}_marked_dup_metrics.txt"
recal_data_table_file="${preprocessing_path}/${sample}_recal_data.table"
recal_bam_file="${bam_path}/${sample}.bam"
recal_bai_file="${bam_path}/${sample}.bai"

# Print out all paths into the slurm log file for easier debugging
echo "-----------------------"
echo "Current working directory: $(pwd)"
echo "Paths for output files:"
echo "fq_R1: $fq_R1"
echo "fq_R2: $fq_R2"
echo "cutadapt_R1: $cutadapt_R1"
echo "cutadapt_R2: $cutadapt_R2"
echo "sam_file: $sam_file"
echo "marked_dup_bam_file: $marked_dup_bam_file"
echo "marked_dup_metrics_file: $marked_dup_metrics_file"
echo "recal_data_table_file: $recal_data_table_file"
echo "recal_bam_file: $recal_bam_file"
echo "recal_bai_file: $recal_bai_file"
echo "-----------------------"
echo ""


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Remove adapter sequences from reads in the raw FASTQ files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load pre-installed modules on O2 required for cutadapt
module load gcc/9.2.0 python/3.8.12 cutadapt/4.1

# create path for cutadapt output if it doesn't exist
mkdir -p $cutadapt_path

# Use cutadapt to remove adapter sequences
if [ ! -f ${adapter_seq_R1} ] || [ ! -f ${adapter_seq_R2} ] ; then
    cutadapt -a ${adapter_seq_R1} -A ${adapter_seq_R2} --minimum-length 20 --cores=20 -o ${cutadapt_R1} -p ${cutadapt_R2} ${fq_R1} ${fq_R2}
fi


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Preprocess sequencing data using GATK's "best practices" for variant discovery
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## load different modules for bwa/gatk
module load gcc/6.2.0 bwa/0.7.15 gatk/4.1.9.0

## create path for output from bwa mem and subsequent steps if it doesn't exist
mkdir -p $preprocessing_path

if [ ! -f ${sam_file} ]; then
    cmd="bwa mem -M -t 20 -R '@RG\tID:"${i}"\tSM:"${sample}"\tPL:${PL}' $reference $cutadapt_R1 $cutadapt_R2 > ${sam_file}"
    echo $cmd
    eval $cmd
fi

if [ -f ${sam_file} ] && [ ! -f ${marked_dup_bam_file} ]; then
    cmd="gatk MarkDuplicatesSpark --input ${sam_file} --output ${marked_dup_bam_file} --tmp-dir ${tmpdir} --reference $reference -M ${marked_dup_metrics_file} --conf 'spark.executor.cores=20'"
    echo $cmd
    eval $cmd
fi

if [ -f ${marked_dup_bam_file} ] && [ ! -f ${recal_data_table_file} ]; then
    cmd="gatk BaseRecalibrator -I ${marked_dup_bam_file} -R $reference --known-sites $polymorphic_sites -O ${recal_data_table_file} --tmp-dir ${tmpdir}"
    echo $cmd
    eval $cmd
fi

if [ -f ${recal_data_table_file} ] && [ ! -f ${recal_bai_file} ]; then
    cmd="gatk ApplyBQSR -I ${marked_dup_bam_file} -R $reference --bqsr-recal-file ${recal_data_table_file} -O ${recal_bam_file} --tmp-dir ${tmpdir}"
    echo $cmd
    eval $cmd
fi


