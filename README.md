# DADA2
# Language: R
# Input: TXT (keyword, value pairs)
# Output: prefix
# Tested with: PluMA 1.1, R 4.0.0
# dada2_1.18.0

PluMA plugin to produce Amplicon Sequence Variants (ASVs) using DADA2 (Callahan et al, 2016)

The plugin accepts as input a TXT file with (keyword, value) pairs.  Three keywords are accepted:

FASTQ: Directory with the input FASTQ sequence files.
DB: Database to use for ASV mapping (i.e. SILVA)
species: Additional species data

Plugin is modeled after the DADA2 tutorial (Callahan, 2020).

The output prefix is used to produce three output files:

<prefix>.fa: Sequence corresponding to each ASV
<prefix>.counts.tsv: Amount of each ASV in each sample
<prefix>.taxonomy.tsv: Taxonomic classification of each ASV
