# Metabarcoding
Metabarcoding Pipeline Slurm/Snakemake tests



  For Crunchomics:
  
  Git clone this repo, activate the environment and run:
  
  ~~~
  git clone https://github.com/PHB-SILS/Metabarcoding
  
  cd Metabarcoding
  
  conda env create --name metabarcoding --file Metaenv.yml
  
  salloc -N 1 -w omics-cn005 --cpus-per-task 30 --mem=30G
  
  conda activate metabarcoding
  
  srun snakemake --cores 30
  ~~~

TO DO:
 - Documentation.
 - Schematic of workflow.
 - Instructions on how, why and when to change trimming parameters.
 - Optimisation (split DADA2 rule and change dada2 script to accommodate more parallelisation).
 - Addition of a couple of R scripts for phyloseq and compositional data corrected phyloseq (diversity measures, PCA etc.) which are nearly ready for integration.
 - DADA2 parameter automation tool (FIGARO or?)
