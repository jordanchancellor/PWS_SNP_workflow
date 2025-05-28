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
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/767/895/GCF_042767895.1_ASM4276789v1/GCF_042767895.1_ASM4276789v1_genomic.fna.gz
     wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/767/895/GCF_042767895.1_ASM4276789v1/GCF_042767895.1_ASM4276789v1_genomic.gtf.gz
     wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/767/895/GCF_042767895.1_ASM4276789v1/GCF_042767895.1_ASM4276789v1_genomic.gff.gz
```
7. Align to downloaded reference genome:
   ```
   sbatch buildindex.sbatch <aligner>
   
   sbatch alignment.sbatch <aligner>
        Aligner options: bwa, bowtie
   ```
8. Get mapping statistics from output files: `sbatch mappingstats.sbatch`

## **Process SAM files for Variant Calling**

9. Clean SAM files based on flags from ValidateSamFile output from `alignment.sbatch`:
    ```
    sbatch cleanbamfiles.sbatch <flags>
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
    ```
    sbatch selectvariants.sbatch <vcf_file> <variant_type> <output_file>
         variant_type options: INDEL, SNP, MIXED, MNP, SYMBOLIC, NO_VARIATION
    ```
22. Check quality measurements on unfiltered SNPs:
```
bcftools query <vcf_file> -f '%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\n' > initial_snps_qualmeasurements.txt
```
23. plot results to view distribution and choose variant filtering thresholds
24. Hard filter variants using either GATK or bcftools based on quality scores above:
```
sbatch filtervariants.sbatch <vcf_file> <filtering_program>
     filter_program options: gatk bcftools
```
25. Repeat steps 22,23 to check new quality metrics on filtered snps
26. Filter for only bi-allelic SNPs, remove monomorphic SNPs, filter on MAF:
```
sbatch filter_biallelicMAF.sbatch <vcf_file> <MAF>
     # MAF is numeric i.e. 0.1, 0.05
```
27. Filter for SNPs only within chromosomes and exclude SNPs within 75bp of beginning/end of chromosomes
- Create BED file from genome fasta:
```
python3 bedfromfasta.py <reference_fasta> <genome_BED_file>
```
- Create BED file of only chromosomes:
```
python3 bedfromfasta.py <reference_fasta> <chromosome_BED_file>
```
- Edit chromosomal BED file to add/subtract 75bp from either end of chromosome:
``` bash
if ($3 - $2 > 150) { 
    start = $2 + 75;    
    end = $3 - 75;      
    print $1"\t"start"\t"end; 
  }
}' <chromosome_BED_file> > <new_chromosome_BED_file>
```
- Index vcf file:
```
bcftools index <vcf_file>
```
- Filter SNPs using bcftools:
```
bcftools view -R ^<new_chromosome_BED_file> -O z -o <new_vcf_file>
```
28. Filter SNPs inside of masked (repetitive) regions
- Mask reference genome using bbmask from BBTools Suite: `sbatch genomemasking.sbatch <reference genome>`
- Create BED file of masked regions (x3: all regions, coding regions, noncoding regions):
```
python3 generate_masked_ranges.py <fasta file | .fa or .fa.gz> > masked_ranges.bed
```
- Index vcf file:
```
bcftools index <vcf_file>
```
- Filter SNPs using bcftools
```
bcftools view -T ^masked_ranges.bed <vcf_file> -O z -o <new_vcf_filename>
```
29. Filter SNPs only in masked, coding regions
- Edit noncoding overlap regions to create bed file of only masked regions which are located in noncoding regions
```
awk -F '\t' '{print $1 "\t" $2 "\t" $3}' masked_in_noncoding_overlap.bed > noncoding_masked_regions.bed
```
- Filter SNPs using bcftools
```
bcftools view -T ^noncoding_masked_regions.bed <vcf_file> -O z -o <new_vcf_file>
```
30. Filter SNPs based on observed heterozygosity
- Calculate heterozygosity per-site: `python3 computeheterozygosity.py <vcf_file>`
- Plot to choose appropriate threshold
- Filter VCF to include only SNPs with heterozygosity greater than specified threshold
```
python3 filtervcf_heterozygosity.py <vcf_file> <heterozygosity_threshold>
```
31.  Filter SNPs based on density across genome: `python3 snp_proximity_filtering.py <vcf_file>`

## **SNP Effect Annotation**

32. Extract annotations of SNPs from reference genome: `python3 getSNPgenomeannotations.py <vcf_file> <reference_genome_gff_file>`
33. Extract flanking sequences of snps: `python3 getSNPflankingregions.py <input_vcf> <referece_genome_fasta_file> <flankingregionupstreamlength> [flankingregiondownstreamlength]`
34. Annotate SNP effects with SNPeff
```
# load SNPeff HPC module
module load snpeff/5e
# copy the snpeff module installation to your new directory so you can edit config files
mkdir snpEff
cd snpEff
scp snpeff/5e ./

# copy .gtf and .fna genome files to data directory
mkdir data
mkdir pvannamei
cd pvannamei
scp <reference_genome_fasta_file> <reference_genome_gtf_file>
gunzip <reference_genome_gtf_file>

# rename files 
mv <reference_genome_fasta_file> sequences.fa
mv <reference_genome_gtf_file> genes.gtf

# add genome to config file
nano snpEff.config
# add this lines plus any other comments you would like to include
# comment out old genomes just to be safe
pvannamei.genome : pvannamei

# create the database 
cd /path/to/snpEff
snpEff build -gtf22 -v pvannamei

# run snpeff
snpEff pvannamei <vcf_file> > <annotated_vcf_filename>
```
### Other scripts included in this repository and not included in the described workflow are used for other SNP analyses. If you have questions about these workflows you can contact me directly 
