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
from lib.IR.get_coverageBed import *
from lib.IR.get_peptide_sequence_RI import *
from lib.IR.select_fasta_candidates import *
from lib.IR.run_netMHC_classI_slurm_part1 import *
from lib.IR.run_netMHC_classI_slurm_part2 import *
from lib.IR.run_netMHCpan_classI_slurm_part1 import *
from lib.IR.run_netMHCpan_classI_slurm_part2 import *

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
        bam_path = "/projects_rg/SCLC_cohorts/George/STAR/George_and_Peifer"
        # bam_path = "/projects_rg/SCLC_cohorts/Rudin/STAR/Rudin_Yokota"
        TPM_threshold = 1
        tumor_specific = True
        transcript_expression_path = "/projects_rg/SCLC_cohorts/George/tables/iso_tpm_George_Peifer_Rudin_Yokota.tab"
        introns_Normal_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/iso_tpm_introns_Rudin_Normal.txt"
        introns_GTEX_path = "/projects_rg/SCLC_cohorts/annotation/chess2.0_assembly_hg19_CrossMap.events_RI_strict.ioe"
        gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.gtf"
        gtf_protein_coding_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        codons_gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.codons.gtf"
        mutations_path = "/projects_rg/babita/TCGA/mutation/mut_pipeline/juanlu_sclc/src_files/SCLC_mutations_sorted.bed.mut.out"
        repeats_path = "/projects_rg/SCLC_cohorts/cis_analysis/tables/hg19_repeats.bed"
        CHESS_SE_path = "/projects_rg/SCLC_cohorts/annotation/chess2.0_assembly_hg19_CrossMap.events_SE_strict.ioe"
        mosea = "/genomics/users/juanluis/Software/MoSEA-master/mosea.py"
        fast_genome = "/genomics/users/juanluis/Software/MoSEA-master/test_files/genome/hg19.fa"
        orfs_scripts = "/genomics/users/juanluis/comprna/MxFinder/extract_orfs.py"
        interpro = "/projects_rg/SCLC_cohorts/soft/interproscan-5.30-69.0/interproscan.sh"
        IUPred = "/projects_rg/SCLC_cohorts/soft/IUPred2A"
        HLAclass_path = "/projects_rg/SCLC_cohorts/tables/PHLAT_summary_ClassI_all_samples.out"
        HLAtypes_path = "/projects_rg/SCLC_cohorts/tables/NetMHC-4.0_HLA_types_accepted.tab"
        HLAtypes_pan_path = "/projects_rg/SCLC_cohorts/tables/NetMHCpan-4.0_HLA_types_accepted.tab"
        netMHC_path = "/users/genomics/juanluis/Software/netMHC-4.0/netMHC"
        netMHC_pan_path = "/users/genomics/juanluis/Software/netMHCpan-4.0/netMHCpan"
        remove_temp_files = True
        coverage_path = "/users/genomics/juanluis/SCLC_cohorts/test2/coverageBed/"
        output_path = "/users/genomics/juanluis/SCLC_cohorts/test2"

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

        # 6. Create the folder, if it doesn't exists
        logger.info("Part6...")
        if not os.path.exists(output_path + "/coverageBed"):
            os.makedirs(output_path + "/coverageBed")
        # Move all the coverage.sorted files to the created directory
        command1="mv "+output_path+"/*coverage_sorted "+output_path + "/coverageBed/"
        os.system(command1)

        # 7. Get the coverage for each exonization
        logger.info("Part7...")
        get_coverageBed(output_path + "/IR_expressed_genes.tab", output_path + "/random_introns.bed",
                        output_path + "/coverageBed", output_path + "/IR_expressed_genes2.tab", True)

        # 8.  Get the peptide sequence associated
        logger.info("Part8...")
        get_peptide_sequence(output_path + "/IR_expressed_genes2.tab", transcript_expression_path, gtf_path, codons_gtf_path,
                             output_path + "/IR_ORF.tab",output_path + "/IR_ORF_sequences.tab", output_path + "/IR_Inerpro.tab",
                             output_path + "/IR_IUPred.tab", mosea, fast_genome, orfs_scripts, interpro,IUPred, remove_temp_files)

        # 9. Filter the significant results
        logger.info("Part9...")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        output_path_aux18 = output_path + "/IR_ORF_filtered.tab"
        output_path_aux19 = output_path + "/IR_ORF_filtered_peptide_change.tab"
        command4="module load R; Rscript "+dir_path+"/lib/IR/filter_results.R "+output_path + "/IR_ORF.tab"+" "+ \
                 output_path + "/IR_ORF_filtered.tab" +" "+output_path + "/IR_ORF_filtered_peptide_change.tab"
        os.system(command4)

        # 10. Select the fasta candidates for being run to the epitope analysis
        logger.info("Part10...")
        #Create the folder, if it doesn't exists
        if not os.path.exists(output_path + "/IR_fasta_files"):
            os.makedirs(output_path + "/IR_fasta_files")
        select_fasta_candidates(output_path + "/IR_ORF_filtered_peptide_change.tab", output_path + "/IR_peptide_sequence.fa", output_path + "/IR_peptide_sequence_filtered.fa", output_path + "/IR_fasta_files")

        #11. Run netMHC-4.0
        logger.info("Part11...")
        if not os.path.exists(output_path + "/IR_NetMHC-4.0_files"):
            os.makedirs(output_path + "/IR_NetMHC-4.0_files")
        run_netMHC_classI_slurm_part1(output_path + "/IR_ORF_filtered_peptide_change.tab", HLAclass_path, HLAtypes_path,
                                      output_path + "/IR_fasta_files",output_path + "/IR_NetMHC-4.0_files", output_path + "/IR_NetMHC-4.0_neoantigens_type_3.tab",
                                      output_path + "/IR_NetMHC-4.0_neoantigens_type_3_all.tab", output_path + "/IR_NetMHC-4.0_neoantigens_type_2.tab",
                                      output_path + "/IR_NetMHC-4.0_neoantigens_type_2_all.tab", netMHC_path)

        #12. Run netMHCpan-4.0
        logger.info("Part12...")
        if not os.path.exists(output_path + "/IR_NetMHCpan-4.0_files"):
            os.makedirs(output_path + "/IR_NetMHCpan-4.0_files")
        run_netMHCpan_classI_slurm_part1(output_path + "/IR_ORF_filtered_peptide_change.tab", HLAclass_path, HLAtypes_pan_path,
                                      output_path + "/IR_fasta_files",output_path + "/IR_NetMHCpan-4.0_files", output_path + "/IR_NetMHCpan-4.0_neoantigens_type_3.tab",
                                      output_path + "/IR_NetMHCpan-4.0_neoantigens_type_3_all.tab", output_path + "/IR_NetMHCpan-4.0_neoantigens_type_2.tab",
                                      output_path + "/IR_NetMHCpan-4.0_neoantigens_type_2_all.tab", netMHC_pan_path)
        logger.info("Wait until all jobs have finished. Then, go on with part3")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()