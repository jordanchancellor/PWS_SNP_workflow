#!/usr/bin/env bash
#SBATCH -J markdupes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=3-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH --array=1-45

# a module-dependent bash script to process alignment files prior to snp calling
# requires that bam files have been cleaned and do not produce any errors when ran through GATK4 ValidateSamFile
# usage: sbatch PWS_SNP_workflow/markdupes.sbatch

# define global variables
THREADS=${SLURM_CPUS_PER_TASK:-12}
wd="/mnt/scratch/users/agr_3216_ifapa/alorenzo/vannamei_snp_workflow"
samples_file="${wd}/samples_file.txt"
sample_id=$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)
outdir="/mnt2/fscratch/users/agr_3216_ifapa/alorenzo/alignment"
reference="${wd}/indices/pvannamei/GCF_042767895.1_ASM4276789v1_genomic.fna" # Genome assembly ASM4276789v1
cleanbam="${outdir}/${sample_id}_sorted_cleaned.bam"
markedsortedbamfile="${outdir}/${sample_id}_sorted_marked.bam"

# Load  environment
module purge
module load samtools/1.21
module load gatk/4.4.0.0

# create output directory (if needed)
if [[ ! -d $outdir ]]; then
	printf "Creating output directory %s\n" "${outdir}"
	mkdir $outdir 
fi

# check GATK4 install
gatk --version
if [[ $? -ne 0 ]]
then
	printf "Error - check GATK4 installation.\n"
	exit 1
fi

# run GATK4's MarkDuplicates
echo "Marking duplicates in ${cleanbam} with GATK's MarkDuplicates."
gatk MarkDuplicates \
-I ${cleanbam} \
-O ${markedsortedbamfile} \
-M ${outdir}/${sample_id}_markdupes.metrics \
--VALIDATION_STRINGENCY LENIENT

if [[ -f ${markedsortedbamfile} ]]
then
	echo "Duplicate marked bam file created for ${sample_id}."
	# run GATK4's ValidateSamFile in summary mode
	echo "Validating ${markedsortedbamfile} with GATK's ValidateSamFile."
	gatk ValidateSamFile --INPUT ${markedsortedbamfile} \
	--REFERENCE_SEQUENCE ${reference} \
	--MODE SUMMARY \
	--OUTPUT ${outdir}/${sample_id}_marked_validate.summary
	# samtools index bam
	echo "Indexing marked, sorted bam file ${markedsortedbamfile}"
	samtools index -@ ${SLURM_CPUS_PER_TASK} \
	-b ${markedsortedbamfile}
else
	echo "Error - duplicate marked bam file not created."
fi
