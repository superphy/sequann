
import os
import dotenv
import glob


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

#FASTAFILES = glob.glob(FASTADIR + '/*.fasta')

rule all:
    input: 
        #DATADIR + "staramr/load_pointfinder.SUCCESS",
        DATADIR + "staramr/load_resfinder.SUCCESS"


rule staramr:
    output: 
        DATADIR + "staramr/pointfinder.tsv",
        DATADIR + "staramr/resfinder.tsv"
    params:
        outdir=DATADIR + "staramr/",
        fastaglob='*.fasta',
        fastadir=FASTADIR
    shell:
        """
        cd {params.fastadir}
        rm -rf {params.outdir}
        staramr search -o {params.outdir} --pointfinder-organism salmonella {params.fastaglob}
        """
        

rule load_resfinder:
    input:
        DATADIR + "staramr/resfinder.tsv"
    output: 
        DATADIR + "staramr/load_resfinder.SUCCESS"
    params:
        analysis="resfinder"
    script:
        "src/loader.py"


# rule load_pointfinder:
#     input:
#         DATADIR + "staramr/pointfinder.tsv"
#     output: 
#         DATADIR + "staramr/load_pointfinder.SUCCESS"
#     params:
#         analysis="pointfinder"
#     script:
#         "src/loader.py"




