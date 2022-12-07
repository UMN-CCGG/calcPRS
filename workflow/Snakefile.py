# load modules
import glob
import os
import subprocess
import pdb
import shutil
import yaml
import pandas as pd
import numpy as np

# Get the date
from datetime import datetime
i = datetime.now()
TIME = i.strftime('%Y-%m-%d')

# specify config file
configfile = "workflow/config.yml"

# parse config.yml file
OUT = config['outname']
MODEL = config['models']
GENOTYPE = config['genotypes']
PYTHON = config['python']
CROSSMAP = config['CrossMap']['exec']

# are the liftover chains provided, or do we need to download -- downloading will take some debugging for now (12/7/2022)
CHAIN_TO_19 = config['CrossMap']['chain_to_hg19']
CHAIN_TO_38 = config['CrossMap']['chain_to_hg38']
if config['CrossMap']['download_chains'] == 'true':
        CHAIN_TO_19 = f"{os.getcwd()}/CrossMapRefData/{os.path.basename(CHAIN_TO_19)}"
        CHAIN_TO_38 = f"{os.getcwd()}/CrossMapRefData/{os.path.basename(CHAIN_TO_38)}"
    else:
        assert (CHAIN_TO_19 != "none" and CHAIN_TO_38 != "none"), "Must provide both liftover chain files in workflow/config.yml file"
        assert os.path.exists(config['CrossMap']['chain_to_hg19']), f"Did not find liftover chain: {config['CrossMap']['chain_to_hg19']}"
        assert os.path.exists(config['CrossMap']['chain_to_hg38']), f"Did not find liftover chain: {config['CrossMap']['chain_to_hg38']}"

# make subdirectories
if not os.path.exists("CrossMapRefData"): os.mkdir("CrossMapRefData")
if not os.path.exists("accessory"): os.mkdir("accessory")
if not os.path.exists("OandE"): os.mkdir("OandE")

for modl in MODEL.keys(): # make a subdirectory for each model
	if not os.path.exists(modl):
		os.mkdir(modl)
		os.mkdir(f"{modl}/tmp")
	for geno in GENOTYPE.keys(): # make a subdirectory for each genotype dataset
		if not os.path.exists(f"{modl}/{geno}"):
			os.mkdir(f"{modl}/{geno}")
			os.mkdir(f"{modl}/{geno}/tmp")

# run locally or submit to cluster depending on config.yml
run_specs = {'__default__':""} # set default to empty string in case of local run
if config['run_settings']['local_run'] == 'true':
	localrules: all # add rule names here once written
else:
	localrules: all
	assert config['run_settings']['scheduler'] in ['pbs', 'slurm'], print(f"\nThe scheduler specified in run_settings -> scheduler: \'{config['run_settings']['scheduler']}\' is invalid.\n")
    assert os.path.exists(config['run_settings']['cluster_config']), print(f"\nMust provide a cluster configuration file under run_settings -> cluster_config in the config.yml file\nFile \'{config['run_settings']['cluster_config']}\' does not exist\n")
    clusterfile = yaml.load(open(config['run_settings']['cluster_config']), Loader=yaml.FullLoader)

# Conditionally set expected outputs of various steps depending on flags in config.yml
def get_all_inputs(wildcards):
	input_list = [f"{OUT}-report.pdf"]
	input_list += expand(f"{{model}}/{{geno}}_final_model", model=MODEL.keys(), geno=GENOTYPE.keys())

rule all: #Main function that determines which files snakemake will attempt to produce upon running
	input: get_all_inputs

rule check_model_for_duplicates: # check model for duplicated sites
	output: "{model}/tmp/duplicated.snp"
	params:
		in_model = lambda wildcards: MODEL[wildcards.model]['model_file']
		prefix = "{model}/tmp/duplicated"
	run:
		shell(f"cat {{params.in_model}} | awk 'FNR>1{ print $1":"$2":"$3":"$4}' | sort | uniq -d  > {{params.prefix}}.snp")

rule make_bed: # make bed file from published PGS model
	output: "{model}/tmp/{model}.bed"
	params:
		in_model = lambda wildcards: MODEL[wildcards.model]['model_file']
		prefix = "{model}/tmp/{model}"
	run:
		shell(f"cat {{params.in_model}} | awk 'FNR>1{print "chr"$1, $2-1, $2}' OFS='\t' > {{params.prefix}}.bed")

rule cross_map: # Convert coordinates, if needed
	input: CHAIN_TO_19, CHAIN_TO_38
	output: f"{{model}}/tmp/{{model}}_lifted.bed"
	params: 
		in_model_build = lambda wildcards: MODEL[wildcards.model]['model_build']
		in_geno_build = lambda wildcards: GENOTYPE[wildcards.geno]['geno_build']
	run:
		if {{in_model_build}} != {{in_geno_build}}:
			if {{in_geno_build}} == "hg38":
				shell(f"{CROSSMAP} bed {CHAIN_TO_38} {{model}}/tmp/{{model}}.bed {{output}}")
			if {{in_geno_build}} == "hg19":
				shell(f"{CROSSMAP} bed {CHAIN_TO_19} {{model}}/tmp/{{model}}.bed {{output}}")
		else: # just cp not lifted bed over to 
			shell(f"cp {{model}}/tmp/{{model}}.bed {{output}}")	