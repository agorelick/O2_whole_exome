# 0. Introduction to the O2 cluster
O2 is a platform for Linux-based high-performance computing at Harvard Medical School. This is the main compute cluster used by the Naxerova Lab at HMS. For reasons of HIPAA compliance and practicality, all “raw” sequencing data should be stored and processed on O2. (Note, “raw” sequencing data refers to sequencing-read-level files like FASTQs, BAMs, SAMs, etc.) Raw sequencing data should never be stored on Dropbox, Google Drive, or other cloud platforms. 

Common operations on raw sequencing data (like mutation calling, copy number profiling, etc.) typically also require large reference data files, which should also be kept on O2. Such reference files may include a version of the human reference genome; and databases of human genetic variation (e.g., common SNPs from gnomAD, dbSNP, 1000 genomes, etc.) As a rule of thumb, all sequencing data that is not processed to the point that you would load into R or open in Excel on a laptop should be worked with on O2. 

To use O2, you will need to register for an account. Instructions are on the O2 wiki: https://harvardmed.atlassian.net/wiki/spaces/O2/overview (See: How to request an O2 account). Requests for new accounts usually take a few days to go through. Note, the wiki for O2 is very well documented and should be your first destination for general information.

## Once you have an O2 account
*	you will need to request that your O2 user be added to the “naxerova” user group. This is necessary for your new user to access the lab’s shared directory on O2 for big data files, which is located at /n/data1/hms/genetics/naxerova/lab/. To request this, submit a ticket on the HMS IT webite (https://harvardmed.service-now.com/stat). _Example ticket description:_ "Hello, I am a postdoc in the Naxerova lab in Department of Genetics at HMS. My account on o2 was recently set up (my username is alg2264). Please add my user to the “naxerova” group, so that I can access our shared data directory at /n/data1/hms/genetics/naxerova/lab. Thank you."
* Create a “scratch” directory (a folder designed for many/large temporarily files that should not be backed up). To create your own scratch directory on O2, please run the following command from a login node (not an interactive job): /n/cluster/bin/scratch_create_directory.sh See: https://harvardmed.atlassian.net/wiki/spaces/O2/pages/2652045313/Scratch+Storage


# 1. QC your FASTQ files 
Before pre-processing your FASTQ files, you can generate a QC report for each file using the commonly used tool FASTQC (www.bioinformatics.babraham.ac.uk/projects/fastqc/). This URL also provides examples of "good" reports and "bad" reports with common problems to aid in your interpretation. To run FASTQC on your data, use the `run_fastqc.sh` script included in this directory. Copy this file to a new file for your current data. Make the following edits:

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

When FASTQC has finished, copy the resulting directory to your local computer to view the results. From a terminal on your local computer, use `scp` to copy the entire fastqc directory. 
```
# run on your local laptop/desktop
# alg2264 is my username on o2, replace with yours
scp -r alg2264@transfer.rc.hms.harvard.edu:/n/data1/hms/genetics/naxerova/lab/alex/example_patient/fastqc example_fastqc
```
Then you can open the fastqc directory (here, `example_fastqc`) and view a file's results (.html files) in a web browser. For information about how to interpret these results, see: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
