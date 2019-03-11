"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu
Neoskipping_ePydoor.py: get significant neoskipping events
"""

import os

from lib.Neoskipping.extract_neoskipping_junctions import *
from lib.Neoskipping.extract_neoskipping_junctions_Intropolis import *
from lib.Neoskipping.check_mutations_nearby import *
from lib.Neoskipping.filter_neoskipping import *
from lib.Neoskipping.filter_neoskipping_CHESS import *

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

        readcounts_path = "/projects_rg/SCLC_cohorts/Hugo/STAR/readCounts.tab"
        transcript_expression_path = "/projects_rg/SCLC_cohorts/George/tables/iso_tpm_George_Peifer_Rudin_Yokota.tab"
        gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        codons_gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.codons.gtf"
        conversion_names = "/projects_rg/SCLC_cohorts/tables/Ensembl_gene_conversion.txt"
        max_length = 500
        threshold = 5
        threshold2 = 10
        repeats_path = "/projects_rg/SCLC_cohorts/cis_analysis/tables/hg19_repeats.bed"
        mutations_path = "/projects_rg/babita/TCGA/mutation/mut_pipeline/juanlu_sclc/src_files/SCLC_mutations_sorted.bed.mut.out"
        CHESS_A5_path = "/projects_rg/SCLC_cohorts/annotation/chess2.0_assembly_hg19_CrossMap.events_A5_strict.ioe"
        CHESS_A3_path = "/projects_rg/SCLC_cohorts/annotation/chess2.0_assembly_hg19_CrossMap.events_A3_strict.ioe"
        tumor_specific = True
        mosea = "/genomics/users/juanluis/Software/MoSEA-master/mosea.py"
        fasta_genome = "/genomics/users/juanluis/Software/MoSEA-master/test_files/genome/hg19.fa"
        orfs_scripts = "/genomics/users/juanluis/comprna/MxFinder/extract_orfs.py"
        interpro = "/projects_rg/SCLC_cohorts/soft/interproscan-5.33-72.0/interproscan.sh"
        IUPred = "/projects_rg/SCLC_cohorts/soft/IUPred2A"
        HLAclass_path = "/projects_rg/SCLC_cohorts/Hugo/Supplementary/HLA_type_Hugo_formatted.tab"
        HLAtypes_path = "/projects_rg/SCLC_cohorts/tables/NetMHC-4.0_HLA_types_accepted.tab"
        HLAtypes_pan_path = "/projects_rg/SCLC_cohorts/tables/NetMHCpan-4.0_HLA_types_accepted.tab"
        netMHC_path = "/projects_rg/SCLC_cohorts/soft/netMHC-4.0/netMHC"
        netMHC_pan_path = "/projects_rg/SCLC_cohorts/soft/netMHCpan-4.0/netMHCpan"
        remove_temp_files = True
        flag_Rudin = False
        output_path = "/users/genomics/juanluis/SCLC_cohorts/Hugo/epydoor/A5_A3"


        # 1. Identify the junctions that could generate an alternative splice site
        logger.info("Part1...")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        output_path_aux = output_path+"/new_Neoskipping_junctions.tab"
        extract_neoskipping_junctions(readcounts_path, gtf_path, threshold, output_path_aux)

        # 2. Get the tumor specific neoskipping events
        if(tumor_specific):

            # Get also the significant neoskipping from Rudin and Intropolis
            output_Rudin_path_aux = output_path + "/new_Neoskipping_junctions_Rudin_normal_reads.tab"
            readCounts_Rudin_path = "/projects_rg/SCLC_cohorts/Rudin/STAR/v1/normal_readCounts.tab"
            extract_neoskipping_junctions(readCounts_Rudin_path, gtf_path, threshold, output_Rudin_path_aux)

            output_Intropolis_path_aux = output_path + "/new_Neoskipping_junctions_Rudin_normal_reads.tab"
            readCounts_Intropolis_path = "/projects_rg/Annotation/Intropolis/intropolis.v1.hg19.filtered.tsv"
            extract_neoskipping_junctions_Intropolis(readcounts_path, readCounts_Intropolis_path, gtf_path, threshold, output_Intropolis_path_aux)

            filter_neoskipping(output_path_aux, output_Rudin_path_aux, output_Intropolis_path_aux, output_path+"/new_Neoskipping_junctions.tab", flag_Rudin)


        logger.info("Wait until all jobs have finished. Then, go on with part2")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()