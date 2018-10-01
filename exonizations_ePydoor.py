"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu
exonizations_ePydoor.py: get significat exonizations
"""

from lib.Exonization.extract_exonized_junctions import *
from lib.Exonization.get_reads_exonizations import *
from lib.Exonization.overlap_with_repeats import *
from lib.Exonization.get_significant_exonizations import *
from lib.Exonization.generate_random_intronic_positions import *

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
        threshold = 5
        n_randomizations = 100
        repeats_path = "/projects_rg/SCLC_cohorts/cis_analysis/tables/hg19_repeats.bed"
        output_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2"

        # 1. Identify the junctions that could generate an exonization
        output_path_aux = output_path+"/new_exonized_junctions.tab"
        extract_exonized_junctions(input_path, gtf_path, max_length, output_path_aux)

        # 2. Given the list with the possible exonizations, get the reads associate to each of them
        output_path_aux2 = output_path+"/new_exonized_junctions_reads.tab"
        get_reads_exonizations(output_path_aux, input_path, output_path_aux2)

        # 3. find the overlap between the nex exonizations and repeatitions (RepeatMasker)
        output_path_aux3 = output_path + "/new_exonized_junctions_reads_repeatitions.tab"
        overlap_with_repeats(output_path_aux2, repeats_path, output_path_aux3)

        # 4. given the table of the exonizations with the reads counts,get those that are over a threshold
        output_path_aux4 = output_path + "/exonizations_by_sample.tab"
        get_significant_exonizations(output_path_aux3, threshold, output_path_aux4)

        # 5. generate a number of random position by exonization
        output_path_aux5 = output_path + "/random_exonizations.gtf"
        output_path_aux6 = output_path + "/random_exonizations.bed"
        generate_random_intronic_positions(output_path_aux4, gtf_path, n_randomizations, output_path_aux5, output_path_aux6)

        logger.info("Done. Exiting program.")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()