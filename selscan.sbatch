#!/usr/bin/env bash
#SBATCH -J selscan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=80G
#SBATCH --time=08:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH --array=2-44

# usuage: sbatch PWS_SNP_workflow/selscan.sbatch
# vcf file must be indexed prior to running: bcftools index <vcf_file>

# set global variables
wd="/mnt2/fscratch/users/agr_3216_ifapa/alorenzo/simulations"
chrom_list="${wd}/chromosomes.txt"
vcf="${wd}/biallelic_qualfiltered_snps.vcf.gz"
chrom=$(cat $chrom_list | sed -n ${SLURM_ARRAY_TASK_ID}p)
shapeit4_path="/mnt2/fscratch/users/agr_3216_ifapa/alorenzo/simulations/shapeit4-4.2.2/bin"
selscan_path="/mnt2/fscratch/users/agr_3216_ifapa/alorenzo/simulations/selscan/src"

# make output directories
mkdir -p ${wd}/phased ${wd}/maps ${wd}/ihs

module purge
module load bcftools/1.16
module load java/jre-8u301

# subset VCF per chromosome
bcftools view -r $chrom $vcf -Oz -o ${wd}/${chrom}.vcf.gz
bcftools index ${wd}/${chrom}.vcf.gz

# generate mock recombination map (1e-8 per bp per gen)
bcftools query -f '%POS\n' ${wd}/${chrom}.vcf.gz > ${wd}/maps/${chrom}_positions.txt
# awk '{pos=$1; print pos/1000, pos*1e-8, pos}' ${wd}/maps/${chrom}_positions.txt > ${wd}/maps/${chrom}.map
awk -v chr=${chrom} '{pos=$1; print chr, "snp"NR, pos*1e-8, pos}' OFS="\t" ${wd}/maps/${chrom}_positions.txt > ${wd}/maps/${chrom}.map
bgzip ${wd}/maps/${chrom}.map

# module purge
# module load boost/1.80
# module load htslib/1.16_exa

# phase subsetted with SHAPEIT4
# ${shapeit4_path}/shapeit4.2 --input ${wd}/${chrom}.vcf.gz \
#          --region ${chrom} \
#          --map ${wd}/maps/${chrom}.map.gz \
#          --output ${wd}/phased/${chrom}.phased.vcf.gz \
#          --log ${wd}/logs/${chrom}.shapeit4.log \
#          --sequencing \
#          --pbwt-depth 2 \
#          --ibd2-count 20000 \
#          --pbwt-modulo 0.0001 \
#          --use-PS 0 \
#          --thread $SLURM_CPUS_PER_TASK

# phase subsetted vcf with beagle
java -Xmx64g -jar ${wd}/beagle.27Feb25.75f.jar \
     gt=${wd}/${chrom}.vcf.gz \
     chrom=${chrom} \
     out=${wd}/phased/${chrom}.phased \
     nthreads=$SLURM_CPUS_PER_TASK
        
# run selscan iHS
${selscan_path}/selscan --ihs \
        --vcf ${wd}/phased/${chrom}.phased.vcf.gz \
        --map ${wd}/maps/${chrom}.map.gz \
        --out ${wd}/ihs/${chrom}_ihs \
        --threads $SLURM_CPUS_PER_TASK

# normalize iHS
${selscan_path}/norm --ihs \
     --bins 20 \
     --files ${wd}/ihs/${chrom}_ihs.ihs.out

echo "Finished $chrom"


