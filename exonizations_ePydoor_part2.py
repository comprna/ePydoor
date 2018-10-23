"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu
exonizations_ePydoor.py: get significat exonizations
"""

import os

from lib.Exonization.extract_exonized_junctions import *
from lib.Exonization.get_reads_exonizations import *
from lib.Exonization.overlap_with_repeats import *
from lib.Exonization.get_significant_exonizations import *
from lib.Exonization.generate_random_intronic_positions import *
from lib.Exonization.get_coverageBed import *
from lib.Exonization.check_mutations_nearby import *

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
        coverage_path = "/projects_rg/SCLC_cohorts/coverageBed/"
        gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        max_length = 500
        threshold = 5
        n_randomizations = 100
        mutations_path = "/projects_rg/babita/TCGA/mutation/mut_pipeline/juanlu_sclc/src_files/SCLC_mutations_sorted.bed.mut.out"
        repeats_path = "/projects_rg/SCLC_cohorts/cis_analysis/tables/hg19_repeats.bed"
        output_path = "/users/genomics/juanluis/SCLC_cohorts/test"

        # 6. Get the coverage for each exonization
        logger.info("Part7...")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        output_path_aux4 = output_path + "/exonizations_by_sample.tab"
        output_path_aux6 = output_path + "/random_exonizations.bed"
        output_path_aux7 = output_path + "/exonizations_by_sample_coverage.tab"
        get_coverageBed(output_path_aux4, output_path_aux6, coverage_path, output_path_aux7)

        # 7. check if in the exonizations there are mutations nearby
        logger.info("Part8...")
        output_path_aux8 = output_path + "/exonizations_by_sample_coverage_mut.tab"
        check_mutations_nearby(output_path_aux7, mutations_path, 200, output_path_aux8)

        # 8. Separate between mutated and non-mutated cases
        logger.info("Part9...")
        output_path_aux9 = output_path + "/mutated_exonizations.tab"
        output_path_aux10 = output_path + "/non_mutated_exonizations.tab"

        command2="Rscript "+dir_path+"/Exonization/separate_mutated_cases.R "+output_path_aux8+" "+output_path_aux9+" "+output_path_aux10
        # print(command2)
        os.system(command2)

        # 9. Join the mutated and non_mutated cases
        logger.info("Part9...")
        output_path_aux11 = output_path + "/all_exonizations.tab"
        command3="cat "+output_path_aux9+" > "+output_path_aux11+";tail -n+2 "+output_path_aux10+" >> "+output_path_aux11
        os.system(command3)

        # 10. Filter the significant results
        logger.info("Part10...")
        output_path_aux12 = output_path + "/all_exonizations_filtered.tab"
        output_path_aux13 = output_path + "/all_exonizations_filtered_peptide_change.tab"
        command4="Rscript "+dir_path+"/Exonization/filter_results.R "+output_path_aux11+" "+output_path_aux12+" "+output_path_aux13
        os.system(command4)


        logger.info("Done. Exiting program.")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()