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
        # bam_path = "/projects_rg/SCLC_cohorts/George/STAR/George_and_Peifer"
        bam_path = "/projects_rg/SCLC_cohorts/Rudin/STAR/Rudin_Yokota"
        coverage_path = "/projects_rg/SCLC_cohorts/coverageBed/"
        gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        max_length = 500
        threshold = 5
        n_randomizations = 100
        mutations_path = "/projects_rg/babita/TCGA/mutation/mut_pipeline/juanlu_sclc/src_files/SCLC_mutations_sorted.bed.mut.out"
        repeats_path = "/projects_rg/SCLC_cohorts/cis_analysis/tables/hg19_repeats.bed"
        output_path = "/users/genomics/juanluis/SCLC_cohorts/test"


        # 1. Identify the junctions that could generate an exonization
        logger.info("Part1...")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        output_path_aux = output_path+"/new_exonized_junctions.tab"
        extract_exonized_junctions(readcounts_path, gtf_path, max_length, output_path_aux)

        # 2. Given the list with the possible exonizations, get the reads associate to each of them
        logger.info("Part2...")
        output_path_aux2 = output_path+"/new_exonized_junctions_reads.tab"
        get_reads_exonizations(output_path_aux, readcounts_path, output_path_aux2)

        # 3. find the overlap between the nex exonizations and repeatitions (RepeatMasker)
        logger.info("Part3...")
        output_path_aux3 = output_path + "/new_exonized_junctions_reads_repeatitions.tab"
        overlap_with_repeats(output_path_aux2, repeats_path, output_path_aux3)

        # 4. given the table of the exonizations with the reads counts,get those that are over a threshold
        logger.info("Part4...")
        output_path_aux4 = output_path + "/exonizations_by_sample.tab"
        get_significant_exonizations(output_path_aux3, threshold, output_path_aux4)

        # 5. generate a number of random position by exonization
        logger.info("Part5...")
        output_path_aux5 = output_path + "/random_exonizations.gtf"
        output_path_aux6 = output_path + "/random_exonizations.bed"
        generate_random_intronic_positions(output_path_aux4, gtf_path, n_randomizations, output_path_aux5, output_path_aux6)

        # 6. Run coverageBed on the samples in the cluster
        logger.info("Part6...")
        command1="for sample in $(ls "+bam_path+"/*/*.sorted.bam | cut -d\"/\" -f7 | cut -d\"_\" -f1 | cut -d\".\" -f1 | sort | uniq );do " \
                "echo \"Processing file $sample: \"$(date); sbatch -J $(echo $sample)_coverageBed "+dir_path+"/coverageBed.sh "+bam_path+"/$(echo $sample)/*.sorted.bam " \
                 " "+output_path_aux6+" "+output_path+"/$(echo $sample).coverage_sorted;done"
        # command1 = "for sample in $(ls " + bam_path + "/*.sorted.bam | cut -d\"/\" -f7 | cut -d\"_\" -f1 | cut -d\".\" -f1 | sort | uniq );do " \
        #                                               "echo \"Processing file $sample: \"$(date); sbatch -J $(echo $sample)_coverageBed " + dir_path + "/coverageBed.sh " + bam_path + "/$(echo $sample).sorted.bam " \
        #                                                " " + output_path_aux6 + " " + output_path + "/$(echo $sample).coverage_sorted;done"
        os.system(command1)
        logger.info("Wait until all jobs have finished. Then, go on with part2")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()