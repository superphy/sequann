
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


#################################################################################
# RULES                                                                         #
#################################################################################

GENOMES, = glob_wildcards(FASTADIR + "/{genome}.fasta")

rule load_resfinder:
    input:  expand(DATADIR + "annotations/resfinder/{genome}_resfinder.txt", genome=GENOMES)
    params:
        inputdir=DATADIR + "annotations/resfinder/",
        analysis="resfinder"
    script:
        "src/loader.py"


rule run_resfinder:
    input:
        FASTADIR + "{genome}.fasta",
    output: DATADIR + "annotations/resfinder/{genome}_resfinder.txt"
    params:
        blastdb=os.environ.get('RESFINDERBLASTDB'),
        analysis="resfinder"
    threads:
        8
    script:
        "src/runner.py"
