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

        # Create unique index
        index1 = pm.IndexModel([("type", pm.DESCENDING), ("genome", pm.DESCENDING), ("contig", pm.DESCENDING),
            ("subject", pm.DESCENDING), ("qstart", pm.DESCENDING), ("qend", pm.DESCENDING)], 
            name="sa_unique_i1", unique=True)

        if not "sa_unique_i1" in self.collection.index_information():
            self.collection.create_indexes([index1])

        if logger:
            self.logger = logger
        else:
            self.logger = logging.getLogger('sequann.src.annot.AnnotDB')



    def size(self):

        return self.collection.count()


    def _valid_document(self, docdict):
       # Check document contains required keywords

        for k in ('type', 'contig', 'genome', 'subject', 'length', 'qstart', 'qend', 
            'qlen', 'sstart', 'send', 'slen'):

            if k not in docdict:
                return 0

        return 1


    def load(self, docdicts):
        """Insert list of documents into mongodb

        Args:
            docdicts(list): List of dictionaries
           
        """

        for d in docdicts:
            if not self._valid_document(d):
                raise ValueError("Invalid format: {}".format(d))

        try:
            result = self.collection.insert_many(docdicts)
            return(result.inserted_ids)
        except pm.errors.BulkWriteError as e:
            self.logger.warn("WriteError: {}".format(e.details))
        
        return(None)


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


def staramr_resfinder_to_json(hit_dict):
    """Convert resfinder hit to MongoDB document

    Args:
        hit_dict(dict): A dictionary with the following keys:
            ('isolate', 'gene', 'pident', 'poverlap', 'len_frac', 'start', 'end', 'contig', 'accession')   
       
    """

    length, slen = [int(s) for s in hit_dict['len_frac'.split('/')]]

    document = {
        'type': 'resfinder',
        'contig': hit_dict['contig'],
        'genome': hit_dict['genome'],
        'subject': "{}|{}".format(hit_dict['gene'],hit_dict['accession']),
        'qstart': int(hit_dict['start']),
        'qend': int(hit_dict['start']),
        'pident': float(hit_dict['pident']),
        'length': length,
        'slen': slen,
        'sstart': None,
        'send': None,
        'overlap': float(hit_dict['poverlap'])
    }

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



