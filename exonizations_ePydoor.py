"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu
MxFinder.py: main script. Check github for more details (https://github.com/JLTrincado/MxFinder)
"""

import sys
import time
import re

from argparse import ArgumentParser, RawTextHelpFormatter

import logging, sys, os
import subprocess
from lib import *

# create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# create console handler and set level to info
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)


def main():
    try:

        logger.info("Starting execution")

        # input_path = sys.argv[1]
        # gtf_path = sys.argv[2]
        # max_length = sys.argv[3]
        # output_path = sys.argv[4]

        input_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering/readCounts_George_Peifer_Rudin_Yokota.tab"
        gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        max_length = 500
        output_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2"

        # 1. Obtain the exon coordinates from the GTF. Format it in a bed format
        output_path_aux = output_path+"/new_exonized_junctions.tab"
        extract_exonized_junctions(input_path, gtf_path, max_length, output_path_aux)

        logger.info("Done. Exiting program.")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()