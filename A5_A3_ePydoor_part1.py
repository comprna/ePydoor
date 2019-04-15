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
from lib.A5_A3.check_mutations_nearby import *
from lib.A5_A3.filter_exonizations import *
from lib.A5_A3.filter_exonizations_CHESS import *
from lib.A5_A3.get_peptide_sequence import *
from lib.A5_A3.select_fasta_candidates import *
from lib.A5_A3.run_netMHC_classI_slurm_part1 import *
from lib.A5_A3.run_netMHCpan_classI_slurm_part1 import *

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

        readcounts_path = "/projects_rg/SCLC_cohorts/Snyder/STAR/readCounts.tab"
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
        interpro = "/soft/EB_repo/bio/sequence/programs/noarch/interproscan/5.33-72.0/interproscan.sh"
        IUPred = "/projects_rg/SCLC_cohorts/soft/IUPred2A"
        HLAclass_path = "/projects_rg/SCLC_cohorts/Snyder/Suppl/HLA_type_Snyder_formatted.tab"
        HLAtypes_path = "/projects_rg/SCLC_cohorts/tables/NetMHC-4.0_HLA_types_accepted.tab"
        HLAtypes_pan_path = "/projects_rg/SCLC_cohorts/tables/NetMHCpan-4.0_HLA_types_accepted.tab"
        netMHC_path = "/projects_rg/SCLC_cohorts/soft/netMHC-4.0/netMHC"
        netMHC_pan_path = "/projects_rg/SCLC_cohorts/soft/netMHCpan-4.0/netMHCpan"
        remove_temp_files = True
        flag_Rudin = False
        output_path = "/users/genomics/juanluis/SCLC_cohorts/Snyder/epydoor/A5_A3"


        # 1. Identify the junctions that could generate an alternative splice site
        logger.info("Part1...")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        output_path_aux = output_path+"/new_A5_A3_junctions.tab"
        extract_exonized_junctions(readcounts_path, gtf_path, max_length, output_path_aux)

        # 2. Given the list with the possible A5_A3, get the reads associate to each of them
        logger.info("Part2...")
        output_path_aux2 = output_path+"/new_A5_A3_junctions_reads.tab"
        get_reads_exonizations(output_path_aux, readcounts_path, output_path_aux2)

        # 3. find the overlap between the nex A5_A3 and repeatitions (RepeatMasker)
        logger.info("Part3...")
        output_path_aux3 = output_path + "/new_A5_A3_junctions_reads_repeatitions.tab"
        overlap_with_repeats(output_path_aux2, repeats_path, output_path_aux3)

        # 4. given the table of the A5_A3 with the reads counts,get those that are over a threshold
        logger.info("Part4...")
        get_significant_exonizations(output_path_aux3, threshold, output_path + "/A5_A3_by_sample.tab")

        # 5. for applying some filtering on the list of A5_A3 junctions, we are gonna compare the readcounts for each
        # junction against other new junctions associated to the same gene
        logger.info("Part5...")
        compare_reads_random_junctions(output_path + "/A5_A3_by_sample.tab", readcounts_path, gtf_path, output_path + "/A5_A3_by_sample_coverage.tab")

        # 6. Check if in the A5_A3 there are mutations nearby
        logger.info("Part6...")
        check_mutations_nearby(output_path + "/A5_A3_by_sample_coverage.tab", mutations_path, 200, output_path + "/A5_A3_by_sample_coverage_mut.tab")

        # 7. Separate the mutated from the non-mutated cases
        logger.info("Part7...")
        command1="module load R; Rscript "+dir_path+"/lib/A5_A3/separate_mutated_cases.R "+ output_path + "/A5_A3_by_sample_coverage_mut.tab " \
                 + output_path + "/A5_A3_mutated.tab " + output_path + "/A5_A3_non_mutated.tab "
        os.system(command1)

        # 8. Get the tumor specific events
        if(tumor_specific):

            # Get also the significant A5_A3 from Rudin and Intropolis
            logger.info("Part8.1...")
            output_Rudin_path_aux2 = output_path + "/new_A5_A3_junctions_Rudin_normal_reads.tab"
            readCounts_Rudin_path = "/projects_rg/SCLC_cohorts/Rudin/STAR/v1/normal_readCounts.tab"
            get_reads_exonizations(output_path+"/new_A5_A3_junctions.tab", readCounts_Rudin_path, output_Rudin_path_aux2)
            output_Rudin_path_aux3 = output_path + "/new_A5_A3_junctions_Rudin_normal_reads_repeatitions.tab"
            overlap_with_repeats(output_Rudin_path_aux2, repeats_path, output_Rudin_path_aux3)
            output_Rudin_path_aux4 = output_path + "/A5_A3_by_sample_Rudin_normal.tab"
            get_significant_exonizations(output_Rudin_path_aux3, threshold2, output_Rudin_path_aux4)

            logger.info("Part8.2...")
            output_Intropolis_path_aux2 = output_path + "/new_A5_A3_junctions_Intropolis_reads.tab"
            get_reads_exonizations(output_path+"/new_A5_A3_junctions.tab", readcounts_path, output_Intropolis_path_aux2)
            output_Intropolis_path_aux3 = output_path + "/new_A5_A3_junctions_Intropolis_reads_repeatitions.tab"
            overlap_with_repeats(output_Intropolis_path_aux2, repeats_path, output_Intropolis_path_aux3)
            output_Intropolis_path_aux4 = output_path + "/A5_A3_by_sample_Intropolis.tab"
            get_significant_exonizations(output_Intropolis_path_aux3, threshold2, output_Intropolis_path_aux4)

            logger.info("Part8.3...")
            output_Rudin_path_aux4 = output_path + "/A5_A3_by_sample_Rudin_normal.tab"
            output_Intropolis_path_aux4 = output_path + "/A5_A3_by_sample_Intropolis.tab"
            output_path_aux11 = output_path + "/A5_A3_non_mutated_filtered.tab"
            filter_exonizations(output_path + "/A5_A3_non_mutated.tab", output_Rudin_path_aux4, output_Intropolis_path_aux4, output_path_aux11, flag_Rudin)
            output_path_aux12 = output_path + "/A5_A3_non_mutated_filtered2.tab"
            filter_exonizations_CHESS(output_path_aux11, CHESS_A5_path, CHESS_A3_path, output_path_aux12)

            # 9. Join the mutated and non_mutated cases
            logger.info("Part8.4...")
            output_path_aux13 = output_path + "/all_A5_A3.tab"
            command3 = "cat " + output_path + "/A5_A3_mutated.tab" + " > " + output_path_aux13 + ";tail -n+2 " + output_path_aux12 + " >> " + output_path_aux13
            os.system(command3)

        else:

            # 9. Join the mutated and non_mutated cases
            logger.info("Part8...")
            output_path_aux13 = output_path + "/all_A5_A3.tab"
            command3 = "cat " + output_path + "/A5_A3_mutated.tab" + " > " + output_path_aux13 + ";tail -n+2 " + output_path + "/A5_A3_non_mutated.tab" + " >> " + output_path_aux13
            os.system(command3)

        # 10. Get the peptide sequence associated
        logger.info("Part9...")
        get_peptide_sequence(output_path_aux13, transcript_expression_path, gtf_path, codons_gtf_path,
                             output_path + "/A5_A3_peptide_sequence.fa", output_path + "/A5_A3_fasta_sequence.fa",
                             output_path + "/A5_A3_ORF.tab", output_path + "/A5_A3_ORF_sequences.tab", output_path + "/A5_A3_Interpro.tab",
                             output_path + "/A5_A3_IUPred.tab", mosea, fasta_genome, orfs_scripts, interpro,IUPred, remove_temp_files)

        # 11. Filter the relevant results
        command4 = "module load R; Rscript " + dir_path + "/lib/A5_A3/filter_results.R " + output_path + "/A5_A3_ORF.tab " \
                   + conversion_names + " " + output_path + "/A5_A3_ORF_filtered.tab " + output_path + "/A5_A3_ORF_filtered_peptide_change.tab"
        os.system(command4)

        # 12. Select the fasta candidates for being run to the epitope analysis
        logger.info("Part10...")
        # Create the folder, if it doesn't exists
        if not os.path.exists(output_path + "/A5_A3_fasta_files"):
            os.makedirs(output_path + "/A5_A3_fasta_files")
        select_fasta_candidates(output_path + "/A5_A3_ORF_filtered_peptide_change.tab",
                                output_path + "/A5_A3_peptide_sequence.fa", output_path + "/A5_A3_peptide_sequence_filtered.fa",
                                output_path + "/A5_A3_fasta_files")

        # 13. Run netMHC-4.0_part1
        logger.info("Part11...")
        if not os.path.exists(output_path + "/A5_A3_NetMHC-4.0_files"):
            os.makedirs(output_path + "/A5_A3_NetMHC-4.0_files")
        run_netMHC_classI_slurm_part1(output_path + "/A5_A3_ORF_filtered_peptide_change.tab", HLAclass_path, HLAtypes_path,
                                      output_path + "/A5_A3_fasta_files",
                                      output_path + "/A5_A3_NetMHC-4.0_files",
                                      output_path + "/A5_A3_NetMHC-4.0_neoantigens_type_3.tab",
                                      output_path + "/A5_A3_NetMHC-4.0_neoantigens_type_3_all.tab",
                                      output_path + "/A5_A3_NetMHC-4.0_neoantigens_type_2.tab",
                                      output_path + "/A5_A3_NetMHC-4.0_neoantigens_type_2_all.tab",
                                      output_path + "/A5_A3_NetMHC-4.0_junctions_ORF_neoantigens.tab",
                                      netMHC_path)

        # 14. Run netMHCpan-4.0_part1
        logger.info("Part12...")
        if not os.path.exists(output_path + "/A5_A3_NetMHCpan-4.0_files"):
            os.makedirs(output_path + "/A5_A3_NetMHCpan-4.0_files")
        run_netMHCpan_classI_slurm_part1(output_path + "/A5_A3_ORF_filtered_peptide_change.tab", HLAclass_path, HLAtypes_pan_path,
                                         output_path + "/A5_A3_fasta_files",
                                         output_path + "/A5_A3_NetMHCpan-4.0_files",
                                         output_path + "/A5_A3_NetMHCpan-4.0_neoantigens_type_3.tab",
                                         output_path + "/A5_A3_NetMHCpan-4.0_neoantigens_type_3_all.tab",
                                         output_path + "/A5_A3_NetMHCpan-4.0_neoantigens_type_2.tab",
                                         output_path + "/A5_A3_NetMHCpan-4.0_neoantigens_type_2_all.tab",
                                         output_path + "/A5_A3_NetMHCpan-4.0_junctions_ORF_neoantigens.tab",
                                         netMHC_pan_path)

        logger.info("Wait until all jobs have finished. Then, go on with part2")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()