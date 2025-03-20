#!/usr/bin/env bash
#SBATCH -J filtervariants
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# a module-dependent bash script to hard-filter variants using GATK4 or bcftools 
# usage: sbatch PWS_SNP_workflow/filtervariants.sbatch <vcf_file> <filtering_program>
# output is a filtered .vcf.gz file
# change values below if you require filters other than default

# define global variables
THREADS=${SLURM_CPUS_PER_TASK:-12}
variant_file=$1
filter_program=$2
outdir=$(dirname "$variant_file")
QD="2.0"
MQ="40.0"
SOR="3.0"
FS="60.0"
MQRankSum="-12.5"
ReadPosRankSum="-8.0"
DP="250.0"

if [ "$#" -ne 2 ]; then
    echo "Error: input files not specified."
    echo "Usage: sbatch PWS_SNP_workflow/filtervariants.sbatch <vcf_file> <filtering_program>"
    echo "filter_program options: gatk, bcftools"
    exit 1
fi

# create output directory (if needed)
if [[ ! -d $outdir ]]; then
	printf "Creating output directory %s\n" "${outdir}"
	mkdir $outdir 
fi

# Load  environment
module purge

# Hard-filter variants depending on specified filtering program
case "$filter_program" in
    gatk)
        echo "Hard-filtering variants in vcf file ${variant_file} using GATK4"
        module load gatk/4.4.0.0
        gatk --java-options "-Xmx80G" VariantFiltration \
            -V ${variant_file} \
            -filter "QD < ${QD}" --filter-name "QD2" \
            -filter "SOR > ${SOR}" --filter-name "SOR3" \
            -filter "FS > ${FS}" --filter-name "FS60" \
            -filter "MQ < ${MQ}" --filter-name "MQ40" \
            -filter "MQRankSum < ${MQRankSum}" --filter-name "MQRankSum-12.5" \
            -filter "ReadPosRankSum < ${ReadPosRankSum}" --filter-name "ReadPosRankSum-8" \
            -filter "DP > ${DP}" --filter-name "DP250" \
            -O ${outdir}/${filter_program}_qualfiltered_snps.vcf.gz
        ;;
    bcftools)
        echo "Hard-filtering variants in vcf file ${variant_file} using bcftools."
        module load bcftools/1.16
        bcftools filter \
        -e "FS>${FS} || SOR>${SOR} || MQ<${MQ} || MQRankSum<${MQRankSum} || QD<${QD} || ReadPosRankSum<${ReadPosRankSum} || INFO/DP>${DP}" \
        -O z -o ${outdir}/${filter_program}_qualfiltered_snps.vcf.gz ${variant_file}
        ;;
    *)
        echo "Error: Invalid filtering program '$filter_program'."
        echo "Usage: sbatch PWS_SNP_workflow/filtervariants.sbatch <vcf_file> <filtering_program>"
        echo "filter_program options: gatk bcftools"
        exit 1
        ;;
esac