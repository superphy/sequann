#!/usr/bin/env python

"""Run external annotation programs

"""

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


import argparse
import logging
import subprocess
from Bio import SeqIO


config = {
    'resfams_hmm_models': '/home/matt/workspace/l_amr/pangenome/annotate/resfams/Resfams.hmm',
    'resfinder_blastdb': '/home/matt/workspace/l_amr/pangenome/annotate/resfinder/resfinder'
}
logger = None


def run_resfinder(options):
    """Run Blast with resfinder sequence DB

    """
    cmd = ['blastn','-query', options.input,
        '-db', config['resfinder_blastdb'], '-out', options.output,
        '-evalue', '0.0001', '-outfmt', 
        "7 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore",
        '-perc_identity', '90']
    logger.info("Running: {}".format(subprocess.list2cmdline(cmd)))

    subprocess.check_output(cmd)


def run_rgi(options):
    """Run RGI

    """
    subprocess.check_output(['rgi', '-i', options.input, '-o', options.output])


def run_resfams(options):
    """Run HMMer on translated sequences

    """

    # Translate
    translation_file = append_it(options.input, 'translated')
    translate_orfs(options.input, translation_file)

    # Run HMMer
    subprocess.check_output(['hmmsearch', '--cpu', '8', '-E', '0.001', '--tblout', options.output, config['resfams_hmm_models'],
        translation_file])


def append_it(filename, id):
    return "{0}_{2}.{1}".format(*filename.rsplit('.', 1) + [id])

   
def translate_orfs(input, output):
    table = 11
    min_pro_len = 70
    with open(output, 'w') as outfh:
        for record in SeqIO.parse(input, "fasta"):
            record_len = len(record)
            for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
                for frame in range(3):

                    length = 3 * ((record_len-frame) // 3) # Multiple of three
                    proteins = nuc[frame:frame+length].translate(table)

                    pos = frame
                    laststop = frame
                    stopcodon = '*'
                    pro = ''

                    for aa in proteins:

                        pos = pos + 3
                        pro = pro + aa

                        if aa == stopcodon:
                            if len(pro) >= min_pro_len:
                                if strand == 1:
                                    b = laststop
                                    e = pos-1
                                else:
                                    b = record_len-laststop-1
                                    e = record_len-(pos-1)-1

                                outfh.write(">{}|strand={};frame={};relstart={};relend={}\n{}\n".format(record.id, 
                                    strand, frame, b, e, pro))

                            pro = ''
                            laststop = pos

                    

if __name__ == "__main__":
    """Run various programs

    """

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('annotate.annot')

    # Parse command-line args
    parser = argparse.ArgumentParser()
    parser.add_argument('program', choices=['rgi','resfams', 'resfinder'],
        help='Annotation program')
    parser.add_argument('input', help='Panseq panGenome.fasta file')
    parser.add_argument('output', help='Output file')
    
    options = parser.parse_args()

    if options.program == 'rgi':
        run_rgi(options)

    elif options.program == 'resfams':
        run_resfams(options)

    elif options.program == 'resfinder':
        run_resfinder(options)



