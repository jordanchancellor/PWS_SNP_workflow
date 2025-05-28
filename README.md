# Workflow for SNP calling in *P. vannamei*

Here is a step-by-step workflow describing how to call SNPs from raw fastqc WGS reads using the scripts in this repository:

## **Read QC & Trimming**

1. Run intial QC on raw reads using FastQC & MultiQC: `sbatch qc.sh`
2. Trim sequences according to your specific requirements using Trimmomatic: `sbatch trimming.sbatch`
3. Re-run FastQC & MultiQC on trimmed reads: `sbatch qc.sh`
4. Repeat process until reads are appropriate quality for mapping
5. Get trimming statistics from output files: `sbatch trimmingstats.sbatch`

## **Reference Genome Mapping**

6. Download reference sequences from NCBI (Genome assembly ASM4276789v1 at the time of analysis):
```wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/767/895/GCF_042767895.1_ASM4276789v1/GCF_042767895.1_ASM4276789v1_genomic.fna.gz
     wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/767/895/GCF_042767895.1_ASM4276789v1/GCF_042767895.1_ASM4276789v1_genomic.gtf.gz
     wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/767/895/GCF_042767895.1_ASM4276789v1/GCF_042767895.1_ASM4276789v1_genomic.gff.gz
```
7. Align to downloaded reference genome:
   ```sbatch alignment.sbatch <aligner>
        Aligner options: bwa, bowtie
   ```
8. Get mapping statistics from output files: `sbatch mappingstats.sbatch`

## **Process SAM files for Variant Calling**

9. Clean SAM files based on flags from ValidateSamFile output from `alignment.sbatch`:
    ```sbatch cleanbamfiles.sbatch <flags>
         flag example: INVALID_TAG_NM
    ```
10. Mark duplicates: `sbatch markdupes.sbatch`
11. Get duplicate statistics from output files: `sbatch duplicate_stats.sbatch`

## **Variant Calling**

12. Call variants per-sample: `sbatch haplotypecaller.sbatch`
13. Split reference genome into intervals for joint genotyping: `sbatch splitintervals.sbatch`
14. Create file list of gVCF files for import GenomicsDBimport:
    ```
         ls ${outdir}/*g.vcf.gz > gvcf.sample.map
    ```
15. Create file to index SplitIntervals:
    ```
         ls ${intervals_dir}/*.intervals > interval_list
    ```
17. Importing gvcf files into split intervals for joint genotyping: `sbatch genomicsDBimport.sbatch`
18. Joint genotype per-interval gVCF files: `sbatch genotypeGVCFs.sbatch`
19. Create list of per-interval VCFs:
    ```
         ls ${outdir}/*.vcf.gz > vcf_files.list
    ```
20. Merge VCF files: `sbatch mergeVCFs.sbatch`

## **Variant Selection & Filtering**

21. 
