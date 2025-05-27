#!/usr/bin/env python

import sys
import msprime
import numpy as np

# population parameter inputs
Ne = 
mu = 
recomb_rate = 1e-8
genome_length = 1.8e9 
snp_positions = np.loadtxt("snp_positions.txt")
num_replicates = 1000
T = np.zeros(num_replicates)

# Simulate ancestry
ts = msprime.sim_ancestry(
    samples=100,  # number of haploid genomes
    population_size=Ne,
    sequence_length=genome_length,
    recombination_rate=recomb_rate, 
)

# Overlay mutations
for pos in snp_positions:
    # Ensure SNP is within simulated region
    if 0 <= pos < sequence_length:
        tables.sites.add_row(position=pos, ancestral_state="0")
        # Add mutation to one of the nodes
        tables.mutations.add_row(site=len(tables.sites)-1, node=0, derived_state="1")

# Create a new tree sequence
ts_with_snps = tables.tree_sequence()

# Save to file 
ts_with_snps.dump("burnin_with_real_snps.trees")
