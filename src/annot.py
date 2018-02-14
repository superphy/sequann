#!/usr/bin/env python

"""Annotate panseq pangenomes

Functions for running AMR gene predictors and cross-referencing them against the panseq pangenome

"""

__author__ = "Matthew Whiteside"
__copyright__ = "Copyright 2015, Public Health Agency of Canada"
__license__ = "APL"
__version__ = "2.0"
__maintainer__ = "Matthew Whiteside"
__email__ = "matthew.whiteside@phac-aspc.gc.ca"


import logging
import os
import pymongo as pm
import re

from dotenv import find_dotenv, load_dotenv

logger = None

class AnnotDB(object):
    """Storage and Retrieval of Panseq Annotations

    Maps annotations on the unfragmented pangenome sequences to the individual 
    pangenome fragments.

    """

    def __init__(self, dburi):
        """Constructor

        Args:
            dburi(str): MongoDB URI
           
           
        """

        client = pm.MongoClient(dburi)
        db = client['sequann']
        self.collection = db['amr']

        self.collection.createIndex( { user: 1, title: 1, Bank: 1 }, {unique:true} )


    def size(self):

        return self.collection.count()


    def _valid_document(self, docdict):
       # Check document contains required keywords

        for k in ('type', 'contig', 'genome', 'subject', 'qstart', 'qend', 'length', 
            'qlen', 'sstart', 'send', 'slen'):

            if k not in docdict:
                return 0

        return 1


def parse_header(header_str):

    m = re.search(r'^(?P<genome>SRR\d+)\.fasta\|(?P<contig>.+)$', str(header_str))
    if not m:
        raise Exception('Invalid header '+str(header_str))

    return (m.group('genome'), m.group('contig'))


def blast_to_json(blast_hit_dict, doctype='blast'):
    """Convert resfinder BLAST hit to MongoDB document

    Args:
        blast_hit_dict(dict): A dictionary with the following keys:
            ('qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend', 'qlen', 'sstart', 
             'send', 'slen', 'evalue', 'bitscore')   
       
    """

    genome, contig = parse_header(blast_hit_dict['qseqid'])

    document = {
        'type': doctype,
        'contig': contig,
        'genome': genome,
        'subject': blast_hit_dict['sseqid'], 
    }

    # Add the following BLAST-specific keywords+data
    for k in ('pident', 'qstart', 'qend', 'length', 'qlen', 'sstart', 
             'send', 'slen', 'evalue', 'bitscore'):

        document[k] = blast_hit_dict[k]


    return document


if __name__ == "__main__":
    """Database status

    """

    load_dotenv(find_dotenv())

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('sequann.src.annot')

    dburi = os.environ.get('DBURI')
  
    ann = AnnotDB(dburi)
    
    logger.info("Number of Documents: {}".format(ann.size()))



