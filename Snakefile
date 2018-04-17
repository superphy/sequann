
import os
import dotenv


#################################################################################
# FUNCTIONS                                                                     #
#################################################################################

def OPJ(*args):
    path = os.path.join(*args)
    return os.path.normpath(path)


#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_NAME = 'sequann'
PROJECT_DIR = OPJ(workflow.basedir)

# dotenv project variables
dotenv_path = OPJ(PROJECT_DIR, ".env")
dotenv.load_dotenv(dotenv_path)

FASTADIR=os.environ.get('SQANFASTADIR')
DATADIR=os.environ.get('SQANDATA')
DBDIR=os.environ.get('DBDIR')


#################################################################################
# RULES                                                                         #
#################################################################################

GENOMES, = glob_wildcards(FASTADIR + "/{genome}.fasta")

rule all:
    input: DBDIR




# rule load_resfinder:
#     input:  expand(DATADIR + "annotations/resfinder/{genome}_resfinder.txt", genome=GENOMES)
#     params:
#         inputdir=DATADIR + "annotations/resfinder/",
#         analysis="resfinder"
#     script:
#         "src/loader.py"


# rule run_resfinder:
#     input:
#         FASTADIR + "{genome}.fasta",
#     output: DATADIR + "annotations/resfinder/{genome}_resfinder.txt"
#     params:
#         blastdb=os.environ.get('RESFINDERBLASTDB'),
#         analysis="resfinder"
#     threads:
#         8
#     script:
#         "src/runner.py"


rule staramr:
    input:
        GENOMES
    output: 
        DATADIR + "annotations/staramr/pointfinder.tsv",
        DATADIR + "annotations/staramr/resfinder.tsv"
    shell:
        "staramr search -o staramr --pointfinder-organism salmonella input"


rule load_resfinder:
    input:
        FASTADIR + "{genome}.fasta"
    output: 
        DATADIR + "annotations/staramr/resfinder.tsv"
    params:
        analysis="resfinder"
    script:
        "src/runner.py"


rule load_pointfinder:
    input:
        FASTADIR + "{genome}.fasta"
    output: 
        DATADIR + "annotations/staramr/pointfinder.tsv"
    params:
        analysis="pointfinder"
    script:
        "src/runner.py"




