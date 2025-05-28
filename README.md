# Workflow for SNP calling in P. vannamei

Here is a step-by-step workflow describing how to call SNPs from raw fastqc WGS reads using the scripts in this repository:

## **Read QC & Trimming**

1. Run intial QC on raw reads using FastQC & MultiQC: `sbatch qc.sh`
2. Trim sequences according to your specific requirements: `sbatch trimming.sbatch`
3. Re-run FastQC & MultiQC on trimmed reads: `sbatch qc.sh`
4. Repeat process until reads are appropriate quality for mapping
5. Get trimming statistics from output files: `sbatch trimmingstats.sbatch`

## **Reference Genome Mapping**

6. Download reference sequences from NCBI (Genome assembly ASM4276789v1 at the time of analysis):
```wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/767/895/GCF_042767895.1_ASM4276789v1/GCF_042767895.1_ASM4276789v1_genomic.fna.gz
     wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/767/895/GCF_042767895.1_ASM4276789v1/GCF_042767895.1_ASM4276789v1_genomic.gtf.gz
     wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/767/895/GCF_042767895.1_ASM4276789v1/GCF_042767895.1_ASM4276789v1_genomic.gff.gz
```
7. Align to downloaded reference genome: `sbatch alignment.sbatch`
8. Get mapping statistics from output files: `sbatch mappingstats.sbatch`

## **Process SAM files for Variant Calling**

9. Clean SAM files based on flags from ValidateSamFile output from `alignment.sbatch`: `sbatch cleanbamfiles.sbatch`
10. Mark duplicates: `sbatch markdupes.sbatch`
11. Mark duplicates: `sbatch markdupes.sbatch`
12. Get duplicate statistics from output files: `sbatch duplicate_stats.sbatch`

## **Variant Calling**

13. Call variants per-sample: `sbatch haplotypecaller.sbatch`
14. Split reference genome into intervals for joint genotyping: `sbatch splitintervals.sbatch`
15. Importing gvcf files into split intervals for joint genotyping: `sbatch genomicsDBimport.sbatch`
16. Joint genotype per-interval gVCF files: `sbatch genotypeGVCFs.sbatch`
17. Merge VCF files: `sbatch mergeVCFs.sbatch`

## **Variant Selection & Filtering**


