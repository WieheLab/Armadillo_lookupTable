# LookupAnalysis


## LookupAnalyais

## Description

Snakemake pipeline used for generating and using armadillo lookup tables. The tables are generated using ARMADiLLO and saved for future use. The python program get_NumbMuts_VDJ.py is used to with the lookup tables to generate the improbable and probable counts. The step by step progress for a run without a lookup table is shown in the diagram below. 

## Installation

This pipeline requires the python 3.8.1 environment with Snakemake installed. Other requirements are  
Bio==1.3.6  
matplotlib==3.5.1  
numpy==1.20.1  
scipy==1.8.0  
xlrd==2.0.1  
and can be found in the requirements.txt file. These can be loaded using pipe as shown:  
python3 -m pip install -r requirements.txt

## Pipeline Usage

snakemake -c1 --configfile config.yaml

snakemake -c50 -s /datacommons/dhvi/scripts/snakemake_pipelines/LookupTable_analysis/Snakefile  --config sample=${WORKING_DIR}${NAME}.human_VDJ.fxnl.fasta  ucaFile=${UCA_FASTA} ucaMarkup=[]

##### command to run snakemake through the slurm processing

snakemake --cluster="sbatch -p dhvi --mem=16G --job-name=${NAME}.snakemake " -c1 -j ${NumberJob} --latency-wait=300 -s /datacommons/dhvi/scripts/snakemake_pipelines/LookupTable_analysis/Snakefile  --config sample=${WORKING_DIR}${NAME}.human_VDJ.fxnl.fasta  ucaFile=${UCA_FASTA} ucaMarkup=[]

##### Partial Clean - cleans out the generated files but not the lookup table or compiled programs

snakemake partialclean -c4 --configfile config.yaml --config basedir=samplefiles --delete-all-output --dry-run

##### Full Clean - cleans out the compiled programs and related lookup table

snakemake fullclean -c4 --configfile config.yaml --config basedir=samplefiles --delete-all-output --dry-run

## Example scripts

snakemake_example.job - job scrip for running DH270, can be inserted at the end of any job script that provides the correct variable

snakemake_cluster.job - job script to submit the snakemake as a slurm job. Each step will spawn its own slurm job

snakemake_individual.job - job script that is a stand alone script for running the pipeline. It can be sumitted and run in the normal way with the first argument being the sample file and the the second being teh ucafile

## DAG

snakemake -c1 --configfile config.yaml --config basedir=samples --forceall --rulegraph | dot -Tpdf > dag.pdf

![Alt text](Lookup_Analysis.dag.png?raw=true "Diagram of lookup Analysis")

## Programs used in analysis

### Compiled Programs

All compiled programs are written in C++ and include makefiles to compile and clean the programs.

**align_all_to_first**

USAGE:
	align_all_to_first -i [fasta file]  -go [gap open] -ge [gap extend] -matrix [scoring matrix file] -retain_ins_to_file [insertion file] -align_type [default=global,semi-global]

**ARMADiLLO**

USAGE:
	 ARMADiLLO [seq file options] -m [S5F mutability file] -s [S5F substitution file] <opt arguments>

**positional_aa_matrix_maker**

USAGE:
	positional_aa_matrix_makerc -i [aligned fasta] -o [frequency table] -transfac [transfac formatted]

**translater**

USAGE:
	translater -i [fasta] -f [1 (default),2,3] -s [flag  to suppress stop codons]

### Python Programs/scripts

**get_NumbMuts_VDJ.py**

  min requirements : python 3.5
  
USAGE:
	python get_NumbMuts.py UCAfile VDJ DirWithtables
  
**plot_selection.py**

  min requirements : Python 3.8

USAGE:
	plot_selection.py [-h] [-d LISTDIRECTORY [LISTDIRECTORY ...]] [-f LISTFILE [LISTFILE ...]] [-o GRAPHNAME] [-v VERSION] [-c CUTOFF] [--observedcutoff OBSERVEDCUTOFF] [--top TOP]

**plot_seq_mutations.py**

  min requirements : Python 3.8.1

USAGE:
	plot_seq_mutations.py [-h] [-u UCAFILE] -s SEQFILES [SEQFILES ...] [-m MARKUPFILE] [-v VERSION]

**process_Freq_Table.py**

  min requirements : python 3.5
  
usage: process_Freq_Table.py [-h] [-u UCA] [-f INFILE] [-c CONVERT] [-l LOGOFILE] [-o OUTFILE] [--writelist WRITELIST]

## Support

Program is supplied as is.

## Contributing

JSMB - 4/14/2022

## License

MIT liscense

