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

21. Select variants of interest:
    ```sbatch selectvariants.sbatch <vcf_file> <variant_type> <output_file>
         variant_type options: INDEL, SNP, MIXED, MNP, SYMBOLIC, NO_VARIATION
    ```
22. Check quality measurements on unfiltered SNPs:
``` bcftools query snps.vcf.gz -f '%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\n' > initial_snps_qualmeasurements.txt
```
23. plot results to view distribution and choose variant filtering thresholds
24. Hard filter variants using either GATK or bcftools based on quality scores above:
``` sbatch filtervariants.sbatch <vcf_file> <filtering_program>
     filter_program options: gatk bcftools
```
25. Repeat steps 22,23 to check new quality metrics on filtered snps
26. Filter for only bi-allelic SNPs, remove monomorphic SNPs, filter on MAF:
``` sbatch filter_biallelicMAF.sbatch <vcf_file> <MAF>
     # MAF is numeric i.e. 0.1, 0.05
```
27. Filter for SNPs only within chromosomes and exclude SNPs within 75bp of beginning/end of chromosomes
```
```
28. 
