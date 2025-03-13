#!/usr/bin/env bash
#SBATCH -J bowtie
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --time=4-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH --array=1-45

# a module-dependent bash script to map fastq files to reference genome using bowtie2 or bwa
# requires you to have first built your desired index using buildindex.sbatch
# usage: sbatch PWS_SNP_workflow/alignment.sbatch <aligner>

# define global variables
THREADS=${SLURM_CPUS_PER_TASK:-12}
aligner=$1
wd="/mnt/scratch/users/agr_3216_ifapa/alorenzo/vannamei_snp_workflow"
samples_file="${wd}/samples_file.txt"
sample_id=$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)
reads="${wd}/reads/trimmed"
outdir="/mnt2/fscratch/users/agr_3216_ifapa/alorenzo/alignment"
read1="${reads}/${sample_id}_F_trimmed_paired.fq.gz"
read2="${reads}/${sample_id}_R_trimmed_paired.fq.gz"
reference="${wd}/indices/pvannamei/GCF_042767895.1_ASM4276789v1_genomic.fna" # Genome assembly ASM4276789v1
index=$(basename "${reference}")
samfile="${outdir}/${sample_id}.sam"
sortedbamfile="${outdir}/${sample_id}_sorted.bam"

# extract read group information from fastq's
ID="$(zcat $read1 | head -n 1 | cut -d':' -f3,4 | tr ':' '.')"
PU="$ID"
LB="$(zcat $read1 | head -n 1 | awk '{print $2}' | cut -f 4 -d ":")"
PL="ILLUMINA"

# Load  environment
module purge
module load samtools/1.21
module load gatk/4.4.0.0

# check if aligner is specified
if [ -z "$1" ]; then
    echo "Error: No aligner specified."
    echo "Usage: PWS_SNP_workflow/alignment.sbatch <aligner>"
    echo "Options: bwa, bowtie"
    exit 1
fi

# create output directory (if needed)
if [[ ! -d $outdir ]]; then
	printf "Creating output directory %s\n" "${outdir}"
	mkdir $outdir 
fi

# align reads as determined by aligner options
case "$aligner" in
    bowtie)
        echo "Running alignment using bowtie2"
        module load bowtie/2.5.1
        # bowtie2 install checkpoint
        bowtie2 --version || { echo "Error: Bowtie2 not installed correctly."; exit 1; }

        echo "Aligning ${read1} and ${read2} to index ${index} using bowtie2." 

        # align reads using bowtie2
        bowtie2 -p $THREADS --no-discordant \
        -x "${index}" \
        --rg-id "$ID" --rg "SM:$sample_id" --rg "LB:$LB" --rg "PL:$PL" \
        -1 $read1 -2 $read2 | \
        # -S $samfile
        samtools view -bh -f 2 --threads $THREADS | \
        samtools sort -m 30G -@ $THREADS -o ${sortedbamfile}
        echo "Bowtie2 mapping completed. Mapping files can be found in ${outdir}."
        ;;
    bwa)
        echo "Running alignment using bwa mem"
        module load bwa/0.7.17
        # BWA install checkpoint
        bwa 2>&1 | head -n 1 || { echo "Error: BWA not installed correctly."; exit 1; }

        echo "Aligning ${read1} and ${read2} to index ${index} using bwa mem."

        # align reads using bwa mem
        bwa mem -t $THREADS \
        -R "@RG\tID:$ID\tSM:$sample_id\tLB:$LB\tPL:$PL" \
        $index $read1 $read2 > "${samfile}"
        echo "BWA mapping completed. Mapping files can be found in ${outdir}."
        ;;
    *)
        echo "Error: Invalid aligner '$aligner'."
        echo "Usage: PWS_SNP_workflow/alignment.sbatch <aligner>"
        echo "Options: bwa, bowtie"
        exit 1
        ;;
esac

# convert sam files to indexed bam files 
# echo "Converting sam files to sorted bam files with samtools."
# samtools view -bh -f 2 --threads $THREADS ${samfile} | samtools sort -m 30G -@ $THREADS -o ${sortedbamfile}

# get mapping statistics from alignment file
echo "Creating mapping rate file with samtools flagstat from ${sortedbamfile}."
samtools flagstat ${sortedbamfile} > ${outdir}/${sample_id}.mappingstats

# run GATK4's ValidateSamFile in summary mode
echo "Validating ${sortedbamfile} with GATK's ValidateSamFile."
gatk ValidateSamFile --INPUT ${sortedbamfile} \
--REFERENCE_SEQUENCE ${reference} \
--MODE SUMMARY \
--OUTPUT ${outdir}/${sample_id}_initial_validate.summary

# run GATK4's CollectAlignmentSummaryMetrics on alignment files
echo "Collecting alignment summary metrics from ${sortedbamfile} with GATK's CollectAlignmentSummaryMetrics."
gatk CollectAlignmentSummaryMetrics \
-R ${reference} \
-I ${sortedbamfile} \
-O ${outdir}/${sample_id}_alignment_summary.metrics

# run GATK4's CollectWgsMetrics on alignment files
echo "Collecting alignment summary metrics from ${sortedbamfile} with GATK's CollectWgsMetrics."
gatk CollectWgsMetrics \
-R ${reference} \
-I ${sortedbamfile} \
-O ${outdir}/${sample_id}_wgs.metrics


