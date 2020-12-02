# Metabarcoding
Metabarcoding Pipeline Slurm/Snakemake tests



  
  Git clone this repo, activate the environment and run:
  
  ~~~
  git clone https://github.com/PHB-SILS/Metabarcoding
  
  cd Metabarcoding
  
  conda env create --name metabarcoding --file Metaenv.yml
  
  conda activate metabarcoding
  
  snakemake --profile slurm/ --jobs 10
  ~~~

TO DO:
    More documentation.
    Schematic of workflow.
    Instructions on how, why and when to change trimming parameters.
    Optimisation (split DADA2 rule and change dada2 script to accommodate more parallelisation).
    Addition of a couple of R scripts for phyloseq and compositional data corrected phyloseq (diversity measures, PCA etc.) which are nearly ready for integration.
    DADA2 parameter automation tool, replacement for FIGARO https://www.biorxiv.org/content/10.1101/610394v1 (biggest and most optional task).
