# Workflow for SNP calling in P. vannamei

Here is a step-by-step workflow describing how to call SNPs from raw fastqc WGS reads using the scripts in this repository:

1. **Read QC & Trimming**

   - Run intial QC on raw reads using FastQC & MultiQC: `sbatch qc.sh`
   - Trim sequences according to your specific requirements: `sbatch trimming.sbatch`
   - Re-run FastQC & MultiQC on trimmed reads: `sbatch qc.sh`
   - Repeat process until reads are appropriate quality for mapping
   - Get trimming statistics from output files: `sbatch trimmingstats.sbatch`

2. **Reference Genome Mapping**

   - Download reference sequences from NCBI (Genome assembly ASM4276789v1 at the time of analysis):
   - ```wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/767/895/GCF_042767895.1_ASM4276789v1/GCF_042767895.1_ASM4276789v1_genomic.fna.gz
     wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/767/895/GCF_042767895.1_ASM4276789v1/GCF_042767895.1_ASM4276789v1_genomic.gtf.gz
     wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/042/767/895/GCF_042767895.1_ASM4276789v1/GCF_042767895.1_ASM4276789v1_genomic.gff.gz```
   - Align to downloaded reference genome: `sbatch alignment.sbatch`
   - Get mapping statistics from output files: `sbatch mappingstats.sbatch`

3. **Process SAM files for Variant Calling**

   - Clean SAM files based on flags from ValidateSamFile output from `alignment.sbatch`: `sbatch cleanbamfiles.sbatch`
   - Mark duplicates: `sbatch markdupes.sbatch`
   - Mark duplicates: `sbatch markdupes.sbatch`
   - Get duplicate statistics from output files: `sbatch duplicate_stats.sbatch`

4. **Variant Calling**

   - Call variants per-sample: `sbatch haplotypecaller.sbatch`
   - Split reference genome into intervals for joint genotyping: `sbatch splitintervals.sbatch`
   - Importing gvcf files into split intervals for joint genotyping: `sbatch genomicsDBimport.sbatch`
   - Joint genotype per-interval gVCF files: `sbatch genotypeGVCFs.sbatch`
   - Merge VCF files: `sbatch mergeVCFs.sbatch`

5. ****


