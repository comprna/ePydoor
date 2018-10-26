"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu
IR_ePydoor.py: get significat exonizations
"""

import os

from lib.IR.extract_significant_IR import *
from lib.IR.IR_associate_gene_ids import *
from lib.IR.filter_IR import *
from lib.IR.filter_IR_CHESS import *
from lib.IR.generate_random_intronic_positions import *

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

        introns_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/iso_tpm_introns_George_Peifer_Rudin_Yokota.txt"
        TPM_threshold = 1
        tumor_specific = True
        introns_Normal_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/iso_tpm_introns_Rudin_Normal.txt"
        introns_GTEX_path = "/projects_rg/SCLC_cohorts/annotation/chess2.0_assembly_hg19_CrossMap.events_RI_strict.ioe"
        gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.gtf"
        gtf_protein_coding_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        output_path = "/users/genomics/juanluis/SCLC_cohorts/test"


        # readcounts_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering/readCounts_George_Peifer_Rudin_Yokota.tab"
        # # bam_path = "/projects_rg/SCLC_cohorts/George/STAR/George_and_Peifer"
        # bam_path = "/projects_rg/SCLC_cohorts/Rudin/STAR/Rudin_Yokota"
        # coverage_path = "/projects_rg/SCLC_cohorts/coverageBed/"
        # gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        # max_length = 500
        # threshold = 5
        # threshold2 = 10
        # n_randomizations = 100
        # mutations_path = "/projects_rg/babita/TCGA/mutation/mut_pipeline/juanlu_sclc/src_files/SCLC_mutations_sorted.bed.mut.out"
        # repeats_path = "/projects_rg/SCLC_cohorts/cis_analysis/tables/hg19_repeats.bed"
        # output_path = "/users/genomics/juanluis/SCLC_cohorts/test"

        # 1. Get the IR expressed
        extract_significant_IR(introns_path, TPM_threshold, output_path + "/IR_expressed.tab")

        # 2. Obtain the gene ids for the introns
        IR_associate_gene_ids(output_path + "/IR_expressed.tab", gtf_path, output_path + "/IR_expressed_genes.tab")

        # 3. Get the IR tumor specific
        if(tumor_specific):

            #Get the significant introns for the set of normal
            extract_significant_IR(introns_Normal_path, TPM_threshold, output_path + "/IR_expressed_Normal.tab")

            #Filter by a set of Normal
            output_path_filtered = output_path + "/IR_expressed_genes_filtered.tab"
            filter_IR(output_path + "/IR_expressed_genes.tab", output_path + "/IR_expressed_Normal.tab", output_path_filtered)

            # Filter by a set of Normal (GTEX)
            output_path_filtered2 = output_path + "/IR_expressed_genes_filtered2.tab"
            filter_IR_CHESS(output_path_filtered, introns_GTEX_path, output_path_filtered2)

        else:
            output_path_filtered2 = output_path + "/IR_expressed_genes.tab"

        # 4. Generate random positions for each intron
        generate_random_intronic_positions(output_path_filtered2, gtf_protein_coding_path, 100, output_path + "/random_introns.gtf",
                                           output_path + "/random_introns.bed")



        # # 1. Identify the junctions that could generate an exonization
        # logger.info("Part1...")
        # dir_path = os.path.dirname(os.path.realpath(__file__))
        # output_path_aux = output_path+"/new_exonized_junctions.tab"
        # extract_exonized_junctions(readcounts_path, gtf_path, max_length, output_path_aux)
        #
        # # 2. Given the list with the possible exonizations, get the reads associate to each of them
        # logger.info("Part2...")
        # output_path_aux2 = output_path+"/new_exonized_junctions_reads.tab"
        # get_reads_exonizations(output_path_aux, readcounts_path, output_path_aux2)
        #
        # # 3. find the overlap between the nex exonizations and repeatitions (RepeatMasker)
        # logger.info("Part3...")
        # output_path_aux3 = output_path + "/new_exonized_junctions_reads_repeatitions.tab"
        # overlap_with_repeats(output_path_aux2, repeats_path, output_path_aux3)
        #
        # # 4. given the table of the exonizations with the reads counts,get those that are over a threshold
        # logger.info("Part4...")
        # output_path_aux4 = output_path + "/exonizations_by_sample.tab"
        # get_significant_exonizations(output_path_aux3, threshold, output_path_aux4)
        #
        # # 5. Get also the significant exonizations from Rudin and Intropolis
        # logger.info("Part5...")
        # output_Rudin_path_aux2 = output_path + "/new_exonized_junctions_Rudin_normal_reads.tab"
        # readCounts_Rudin_path = "/projects_rg/SCLC_cohorts/Rudin/STAR/v1/normal_readCounts.tab"
        # get_reads_exonizations(output_path_aux, readCounts_Rudin_path, output_Rudin_path_aux2)
        # output_Rudin_path_aux3 = output_path + "/new_exonized_junctions_Rudin_normal_reads_repeatitions.tab"
        # overlap_with_repeats(output_Rudin_path_aux2, repeats_path, output_Rudin_path_aux3)
        # output_Rudin_path_aux4 = output_path + "/exonizations_by_sample_Rudin_normal.tab"
        # get_significant_exonizations(output_Rudin_path_aux3, threshold2, output_Rudin_path_aux4)
        #
        # output_Intropolis_path_aux2 = output_path + "/new_exonized_junctions_Intropolis_reads.tab"
        # get_reads_exonizations(output_path_aux, readcounts_path, output_Intropolis_path_aux2)
        # output_Intropolis_path_aux3 = output_path + "/new_exonized_junctions_Intropolis_reads_repeatitions.tab"
        # overlap_with_repeats(output_Intropolis_path_aux2, repeats_path, output_Intropolis_path_aux3)
        # output_Intropolis_path_aux4 = output_path + "/exonizations_by_sample_Intropolis.tab"
        # get_significant_exonizations(output_Intropolis_path_aux3, threshold2, output_Intropolis_path_aux4)
        #
        # # 6. generate a number of random position by exonization
        # logger.info("Part6...")
        # output_path_aux5 = output_path + "/random_exonizations.gtf"
        # output_path_aux6 = output_path + "/random_exonizations.bed"
        # generate_random_intronic_positions(output_path_aux4, gtf_path, n_randomizations, output_path_aux5, output_path_aux6)
        #
        # # 7. Run coverageBed on the samples in the cluster
        # logger.info("Part7...")
        # command1="for sample in $(ls "+bam_path+"/*/*.sorted.bam | cut -d\"/\" -f7 | cut -d\"_\" -f1 | cut -d\".\" -f1 | sort | uniq );do " \
        #         "echo \"Processing file $sample: \"$(date); sbatch -J $(echo $sample)_coverageBed "+dir_path+"/coverageBed.sh "+bam_path+"/$(echo $sample)/*.sorted.bam " \
        #          " "+output_path_aux6+" "+output_path+"/$(echo $sample).coverage_sorted;done"
        # # command1 = "for sample in $(ls " + bam_path + "/*.sorted.bam | cut -d\"/\" -f7 | cut -d\"_\" -f1 | cut -d\".\" -f1 | sort | uniq );do " \
        # #                                               "echo \"Processing file $sample: \"$(date); sbatch -J $(echo $sample)_coverageBed " + dir_path + "/coverageBed.sh " + bam_path + "/$(echo $sample).sorted.bam " \
        # #                                                " " + output_path_aux6 + " " + output_path + "/$(echo $sample).coverage_sorted;done"
        # os.system(command1)
        # logger.info("Wait until all jobs have finished. Then, go on with part2")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()