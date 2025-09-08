# Fine-mapping pipeline for use with slurm cluster and google cloud 

This repository includes code for finemapping GWAS loci with methods SuSiE and CARMA. 

Generally, *.sh files are used to submit Rscripts to cluster.

## Fine-mapping Steps 

### Step 1

Subset GWAS summary statistics to generate loci-specific files, for given window about lead SNP.

### Step 2

Generate LD matrices for loci from UKBB hail block matrices (requires connecting to cloud).

### Step 3

Back on cluster, run SuSiE and CARMA on loci. 

## Cite

TK

## Maintainer

Sophie Gunn @ sgunn@nygenome.org
TJ Singh lab @ singhlab@nygenome.org

## Acknowledgements

## Release Notes