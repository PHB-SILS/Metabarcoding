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

