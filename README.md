# Chromosomal Rearrangements: Tempo and Mode of Karyotype Evolution in Scarabaeoidea

## Authors
* Sean Chien (schien@tamu.edu)
* Heath Blackmon (coleoguy@gmail.com)

## Summary
This repository contains the data and custom scripts used to analyze chromosomal evolution and rearrangements across the Scarabaeoidea superfamily.

## Software Requirements
* R

## Project Layout & Workflow
To replicate the results, please execute scripts in the following order:

### /scripts
Contains all analysis code.
* `basic.model.R`:

Part 1: Models chromosome evolution under XO and XY systems. Output: `results/chrom_number_model_result.rds`

Part 2: Models chromosome evolution under XO, XY, and NeoXY systems. Output:`results/simple_model_scs.rds`

* `subtree_simple_model.R`:

Family-specific subtree analyses for Scarabaeidae, Passalidae, and Lucanidae. Output: `/results/sub_sca.csv`, `/results/sub_pas.csv`, and `/results/sub_luc.csv`


* `SA_fusion_prob.R`:

SA-fusion analysis. Output: `results/simmap/simmap_*.rds'`, `results/simmap/exp.csv` and `results/simmap/obs.csv`
Generates the SA-fusion representation plots used in the manuscript using Base R (Figure 5).

* `species_level_analyses.R`:

A sensitivity analysis that excludes genus-level taxa to ensure results are robust at the species level. 

Subtree Output: `/results/sub_species_level_sca.csv`, `/results/sub_species_level_pas.csv`, and `/results/sub_species_level_luc.csv`

SA-fusion analysis Output: `results/simmap_species_level/simmap_*.rds'`

Generates the species level only posterior distribution of rates of fusions and fission in three Scarabaeoidae families (Figure S3) and SA-fusion plot used in the manuscript using Base R (Figure S4).

### /data
* `final100trees`: A multiphylo file containing 100 posterior probability trees generated from [BEAST2].
* `SpeciesChromList.csv`: Chromosome data sourced from [Karyotype Databases](https://coleoguy.github.io/karyotypes/index.html).

### /results
This folder contains the primary output files from the analyses. All results can be fully regenerated using the code provided in the /scripts directory.

* `TII_Mesquite.csv`: Results of the Taxonomic Instability Index (TII) calculated via Mesquite.


### /figure
Contains R scripts (Base R) used to generate the figures presented in the manuscript.

`Figure_2A_scs_autosome_dist.R`: Figure 2A. This figure visualizes all 478 available karyotype records from the database, showing the distribution of haploid autosome numbers and sex chromosome system (SCS). 

`Figure_2B_XY-NeoXY.R`: Figure 2B. There are 7 genera that possess both XY and neo-XY SCSs in available karyotype records from the database.

`Figure_3_phylo-chrom-heatmap.R`: Figure 3. Haploid autosome number and species count across genera in Scarabaeidae, Passalidae, and Lucanidae. 

`Figure_4_fusion-fission_rate.R`: Figure 4. Posterior distribution of rates of fusions and fission in three Scarabaeoidae families.

`Figure_S1_TII_plot.R`: Figure S1. Taxonomic instability indices. 

`Figure_S2_all_fusion-fission.R`: Figure S2. Rate of fusion and fission in Scarabaeoidae.
