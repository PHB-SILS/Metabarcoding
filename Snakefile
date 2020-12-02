#########################################
# Snakemake pipeline for Metabarcoding
#########################################


###########
# Libraries
###########
import pandas as pd
#from pytools.persistent_dict import PersistentDict


###############
# Configuration
###############
configfile: "config.yaml" # where to find parameters
WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]

#storage = PersistentDict("mystorage")

units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

# create lists containing the sample names and conditions
SAMPLES = units.index.get_level_values('sample').unique().tolist()
samples = pd.read_csv(config["units"], dtype=str,index_col=0,sep="\t")
ADAPTERFILE = config["fastp"]["adapters"]
samplefile = config["units"]

#get URLs for taxa db (for DADA2 rule)
taxa_train_url = config["refs"]["taxa_train"]
taxa_spec_url = config["refs"]["taxa_spec"]



##From Koes RNAseq pipeline##
# Input functions for rules
#############################

def sample_is_single_end(sample):
    """This function detect missing value in the column 2 of the units.tsv"""
    if "fq2" not in samples.columns:
        return True
    else:
        return pd.isnull(samples.loc[(sample), "fq2"])

def get_fastq(wildcards):
    """ This function checks if the sample has paired end or single end reads
    and returns 1 or 2 names of the fastq files """
    if sample_is_single_end(wildcards.sample):
        return "16S/" + samples.loc[(wildcards.sample), ["fq1"]].dropna()
    else:
        return "16S/" + samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

#def get_trimmed(wildcards):
#    """ This function checks if sample is paired end or single end
#    and returns 1 or 2 names of the trimmed fastq files """
#    if sample_is_single_end(wildcards.sample):
#        return WORKING_DIR + "trimmed/" + samples.loc[(wildcards.sample)] + "_R1_trimmed.fq.gz"
#    else:
#        return WORKING_DIR + "trimmed/" + samples.loc[(wildcards.sample)] + "_R1_trimmed.fq.gz", WORKING_DIR + "trimmed/" + samples.loc[(wildcards.sample)] + "_R2_trimmed.fq.gz"


#################
# Desired outputs
#################
rule all:
    input:
        expand(WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz", sample = SAMPLES),
        #RESULT_DIR + "trimParameters.json"
        RESULT_DIR + "asv.txt",
        RESULT_DIR + "taxa.txt"
    message:
        "Job finished"

#######
# Rules
#######

#####################
# Download Taxa db 
#####################

rule get_taxa_train:
    output:
        "silva_nr99_v138_train_set.fa.gz"
    message:
        "downloading silva train set"
    shell:
        "wget {taxa_train_url}"

rule get_taxa_spec:
    output:
        "silva_species_assignment_v138.fa.gz"
    message:
        "downloading the silva species assignment data"
    shell:
        "wget {taxa_spec_url}"



##################################
# Fastp
##################################

rule fastp:
    input:
        get_fastq,
        ADAPTERFILE
    output:
        fq1  = WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz",
        fq2  = WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz",
        html = RESULT_DIR + "fastp/{sample}.html"
    message:"trimming {wildcards.sample} reads"
    threads: 8
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        sampleName = "{sample}",
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"],
        cut_front_window_size = config["fastp"]["cut_front_window_size"],
        cut_front_mean_quality = config["fastp"]["cut_front_mean_quality"]
    shell:
        "fastp --thread {threads}  --html {output.html} \
        -Q \
        --cut_front \
        --cut_front_window_size {params.cut_front_window_size} \
        --cut_front_mean_quality {params.cut_front_mean_quality} \
        --adapter_fasta {ADAPTERFILE} \
        --in1 {input[0]} --in2 {input[1]} --out1 {output.fq1} --out2 {output.fq2}; \
        2> {log}"

#   --detect_adapter_for_pe \     
##############
#FIGARO
##############

#git clone https://github.com/Zymo-Research/figaro.git
#cd figaro
#pip3 install -r requirements.txt
#python3 figaro.py -i /path/to/fastq/directory -o /path/to/output/files -a [amplicon length] \
#    -f [forward primer length] -r [reverse primer length]


###FIGARO does not work, cannot debug - needs replacing (currently no alternative tools for automating this stage)
#rule FIGARO:
#    input:
#        directory("~/MetabarcodingCOPY/Metabarcoding/16S")
#    output:
#        dir = directory(RESULT_DIR + "FIGARO/"),
#        tp = RESULT_DIR + "FIGARO/trimParameters.json"
#    message: "Optimising DADA2 Parameters"
#    threads: 8
#    shell:
#        "python3 ~/MetabarcodingCOPY/Metabarcoding/figaro/figaro.py -i {input} -o {output.dir} -a 500 -f 39 -r 42 -F eurofins"



#python3 figaro.py -i ~/MetabarcodingCOPY/Metabarcoding/ITS -o ~/MetabarcodingCOPY/Metabarcoding/ITS/ -a 500 -f 39 -r 42 -F PreTrimmed

#python3 ~/MetabarcodingCOPY/Metabarcoding/figaro/figaro.py -i ~/MetabarcodingCOPY/Metabarcoding/16S/ -o ~/MetabarcodingCOPY/Metabarcoding/results/FIGARO/ -a 450 -f 39 -r 42 -F eurofins




#how to get input to work on this rule? need trimmed fastq but 
#Not all output, log and benchmark files of rule contain the same wildcards...
###############
#DADA2
###############

rule DADA2:
    input:
        forward_read = expand(WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz", sample = SAMPLES),
        reverse_read = expand(WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz", sample = SAMPLES),
        TaxAssign = "silva_nr99_v138_train_set.fa.gz",
        SpecAdd = "silva_species_assignment_v138.fa.gz"
    output:
        #filt_fastq_F = WORKING_DIR + "trimmed/filtered/" + "{sample}_F_filt.fastq.gz",
        #filt_fastq_R = WORKING_DIR + "trimmed/filtered/" + "{sample}_R_filt.fastq.gz",
        taxa = RESULT_DIR + "taxa.txt",
        asv = RESULT_DIR + "asv.txt"
    message: "Creating ASV table with DADA2, this will take a while."
    threads: 8
#    script: "scripts/DADAsnake.R"
    shell: "Rscript --vanilla scripts/DADAsnake.R {input.TaxAssign} {input.SpecAdd} {output.taxa} {output.asv}"



################
#Normal phyloseq
################

#rule phylo:
#    input:
#        
#        
#    output:



#############################################
#CoDA corrected phyloseq + diversity measures
#############################################

#rule CODA:
#    input:
        


