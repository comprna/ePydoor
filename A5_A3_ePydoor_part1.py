"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu
A5_A3_ePydoor.py: get significant alternative splice site events
"""

import os

from lib.A5_A3.extract_exonized_junctions import *
from lib.A5_A3.get_reads_exonizations import *
from lib.A5_A3.overlap_with_repeats import *
from lib.A5_A3.get_significant_exonizations import *
from lib.A5_A3.compare_reads_random_junctions import *

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

        # readcounts_path = sys.argv[1]
        # gtf_path = sys.argv[2]
        # max_length = sys.argv[3]
        # output_path = sys.argv[4]

        readcounts_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering/readCounts_George_Peifer_Rudin_Yokota.tab"
        bam_path = "/projects_rg/SCLC_cohorts/George/STAR/George_and_Peifer"
        gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        max_length = 500
        threshold = 5
        threshold2 = 10
        n_randomizations = 100
        repeats_path = "/projects_rg/SCLC_cohorts/cis_analysis/tables/hg19_repeats.bed"
        output_path = "/users/genomics/juanluis/test_Junckey_v2"

        # 1. Identify the junctions that could generate an alternative splice site
        logger.info("Part1...")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        output_path_aux = output_path+"/new_A5_A3_junctions.tab"
        extract_exonized_junctions(readcounts_path, gtf_path, max_length, output_path_aux)

        # 2. Given the list with the possible exonizations, get the reads associate to each of them
        logger.info("Part2...")
        output_path_aux2 = output_path+"/new_A5_A3_junctions_reads.tab"
        get_reads_exonizations(output_path_aux, readcounts_path, output_path_aux2)

        # 3. find the overlap between the nex exonizations and repeatitions (RepeatMasker)
        logger.info("Part3...")
        output_path_aux3 = output_path + "/new_A5_A3_junctions_reads_repeatitions.tab"
        overlap_with_repeats(output_path_aux2, repeats_path, output_path_aux3)

        # 4. given the table of the exonizations with the reads counts,get those that are over a threshold
        logger.info("Part4...")
        output_path_aux4 = output_path + "/A5_A3_by_sample.tab"
        get_significant_exonizations(output_path_aux3, threshold, output_path_aux4)

        # 5. for applying some filtering on the list of A5_A3 junctions, we are gonna compare the readcounts for each
        # junction against other new junctions associated to the same gene
        logger.info("Part5...")
        output_path_aux4 = output_path + "/A5_A3_by_sample_coverage.tab"
        compare_reads_random_junctions(output_path + "/A5_A3_by_sample.tab", readcounts_path, gtf_path, output_path_aux4)

        logger.info("Wait until all jobs have finished. Then, go on with part2")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()