#!/usr/bin/env Rscript

# Usage: Rscript PWS_SNP_workflow/LD.R <input_prefix> <output_dir> <chromosome_list>
# Note: Must have modules R and plink loaded into environment prior to running
# Note: Requires input list of chromosome names

# Load necessary libraries
library(ggplot2)
library(data.table)

# Parse command-line arguments
cmd_args <- commandArgs(trailingOnly = TRUE)
print(cmd_args)

if (length(cmd_args) < 3) {
  stop("Usage: Rscript LD.R <input_prefix> <output_dir> <chromosome_list>")
}

# Define variables from inputs
input_prefix <- cmd_args[1]
output_dir <- cmd_args[2]
chromosome_list_file <- cmd_args[3]

# Read chromosome names from file
if (!file.exists(chromosome_list_file)) {
  stop("Chromosome list file not found:", chromosome_list_file)
}
chromosomes <- readLines(chromosome_list_file)

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to calculate LD using PLINK and plot results
calculate_ld <- function(input_prefix, output_dir, chromosome) {
  # Define PLINK output file name
  plink_out <- file.path(output_dir, paste0("chr_", chromosome, "_ld"))

  # Run PLINK command to calculate LD
  plink_cmd <- paste(
    "plink --bfile", input_prefix,
    "--chr-set 44",
    "--allow-extra-chr",
    "--chr", chromosome,
    "--maf 0.01",
    "--mind 0.02",
    "--hwe 1e-6",
    "--r2",
    "--ld-window-kb 1000",
    "--ld-window 1000",
    "--ld-window-r2 0.1",
    "--out", plink_out
  )
  system(plink_cmd)

  # Read PLINK output file
  ld_file <- paste0(plink_out, ".ld")
  if (!file.exists(ld_file)) {
    message("PLINK output file not found for chromosome ", chromosome)
    return(NULL)
  }
  
  ld_data <- fread(ld_file)

  # Ensure required columns exist
  if (!all(c("BP_A", "BP_B", "R2") %in% colnames(ld_data))) {
    message("Missing required columns in PLINK output for chromosome ", chromosome)
    return(NULL)
  }

  # Calculate distance (kb)
  ld_data[, Distance_kb := abs(BP_B - BP_A) / 1000]

  # Plot LD decay
  plot_file <- paste0(plink_out, "_ld_decay.png")
  p <- ggplot(ld_data, aes(x = Distance_kb, y = R2)) +
    geom_point(alpha = 0.5, color = "blue") +
    geom_smooth(method="loess", formula="y ~ x") +
    labs(title = paste("LD Decay - Chromosome", chromosome),
         x = "Distance (kb)",
         y = expression(r^2)) +
    theme_minimal()

  # Save plot
  ggsave(plot_file, p, width = 8, height = 6)
  message("LD plot saved:", plot_file)
}

# Process each chromosome
for (chr in chromosomes) {
  message("Processing chromosome:", chr)
  calculate_ld(input_prefix, output_dir, chr)
}