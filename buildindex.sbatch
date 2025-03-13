#!/usr/bin/env bash
#SBATCH -J indexbuild
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# a module-dependent bash script build your desired index for reference genome mapping
# usage: sbatch PWS_SNP_workflow/buildindex.sbatch <aligner> <reference>

# check if aligner is specified
if [ "$#" -lt 2 ]; then
    echo "Error: No aligner specified."
    echo "Usage: PWS_SNP_workflow/buildindex.sbatch <aligner> <reference>"
    echo "Options: bwa, bowtie"
    exit 1
fi

# define global variables
THREADS=${SLURM_CPUS_PER_TASK:-12}
aligner=$1
reference=$2

# unzip reference (if necessary) and update variable
if [[ "$reference" == *.gz ]]; then
    echo "Reference genome is compressed. Unzipping..."
    gunzip -c "$reference" > "${reference%.gz}"
    reference="${reference%.gz}"
    echo "Unzipped reference: $reference"
fi

# # define index prefix
index=$(basename "${reference}")

# align reads as determined by aligner options
case "$aligner" in
    bowtie)
        echo "Building Bowtie2 index..."
        module load bowtie/2.5.1
        # bowtie2 install checkpoint
        bowtie2 --version || { echo "Error: Bowtie2 not installed correctly."; exit 1; }

        # check for existence of index files, if none exist, create index
        if [[ ! -f ${index}.1.bt2 ]]; then
            echo "Bowtie2 index not found. Creating index..."
            bowtie2-build --threads $THREADS -f \
            ${reference} ${index}
        fi
        echo "Bowtie index with prefix ${index} built and can be found in ${wd}."
        ;;
    bwa)
        echo "Building BWA index..."
        module load bwa/0.7.17
        # BWA install checkpoint
        bwa 2>&1 | head -n 1 || { echo "Error: BWA not installed correctly."; exit 1; }

        # check for existence of index files, if none exist, create index
        if [[ ! -f ${index}.amb  ]]; then
        bwa index -p ${index} \
        -a bwtsw \
        ${reference}
        fi
        echo "BWA index with prefix ${index} built and can be found in ${wd}."
        ;;
    *)
        echo "Error: Invalid aligner '$aligner'."
        echo "Usage: PWS_SNP_workflow/buildindex.sbatch <aligner> <reference>"
        echo "Options: bwa, bowtie"
        exit 1
        ;;
esac

# Build samtools FASTA index
echo "Building samtools reference index"
module purge
module load samtools/1.21
samtools faidx ${reference}

# Build GATK FASTA dictionary
echo "Building GATK reference dictionary"
module purge
module load gatk/4.4.0.0
gatk CreateSequenceDictionary -R ${reference}
