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
from lib.IR.get_coverageBed_adapter import *
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

        introns_path = "/homes/users/jtrincado/scratch/test_Junckey/iso_tpm_introns_George_Peifer_Rudin_Yokota.txt"
        bam_path = "/homes/users/jtrincado/scratch/test_Junckey/George_and_Peifer"
        # bam_path = "/projects_rg/SCLC_cohorts/Rudin/STAR/Rudin_Yokota"
        TPM_threshold = 1
        tumor_specific = True
        transcript_expression_path = "/homes/users/jtrincado/scratch/test_Junckey/iso_tpm_George_Peifer_Rudin_Yokota.tab"
        introns_Normal_path = "/projects_rg/SCLC_cohorts/cis_analysis/v5/SCLC_v5/tables/iso_tpm_introns_Rudin_Normal.txt"
        introns_GTEX_path = "/projects_rg/SCLC_cohorts/annotation/chess2.0_assembly_hg19_CrossMap.events_RI_strict.ioe"
        gtf_path = "/homes/users/jtrincado/scratch/test_Junckey/Homo_sapiens.GRCh37.75.formatted.gtf"
        gtf_protein_coding_path = "/homes/users/jtrincado/scratch/test_Junckey/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        codons_gtf_path = "/homes/users/jtrincado/scratch/test_Junckey/Homo_sapiens.GRCh37.75.codons.gtf"
        mutations_path = "/projects_rg/babita/TCGA/mutation/mut_pipeline/juanlu_sclc/src_files/SCLC_mutations_sorted.bed.mut.out"
        repeats_path = "/homes/users/jtrincado/scratch/test_Junckey/hg19_repeats.bed"
        CHESS_SE_path = "/homes/users/jtrincado/scratch/test_Junckey/chess2.0_assembly_hg19_CrossMap.events_SE_strict.ioe"
        mosea = "/homes/users/jtrincado/scratch/Software/MoSEA-master/mosea.py"
        fast_genome = "/homes/users/jtrincado/scratch/Software/MoSEA-master/test_files/genome/hg19.fa"
        orfs_scripts = "/homes/users/jtrincado/scratch/Software/MxFinder/extract_orfs.py"
        interpro = "/homes/users/jtrincado/scratch/Software/interproscan-5.30-69.0/interproscan.sh"
        IUPred = "/homes/users/jtrincado/scratch/Software/IUPred2A"
        HLAclass_path = "/homes/users/jtrincado/scratch/test_Junckey/PHLAT_summary_ClassI_all_samples.out"
        HLAtypes_path = "/homes/users/jtrincado/scratch/test_Junckey/NetMHC-4.0_HLA_types_accepted.tab"
        HLAtypes_pan_path = "/homes/users/jtrincado/scratch/test_Junckey/NetMHCpan-4.0_HLA_types_accepted.tab"
        netMHC_path = "/homes/users/jtrincado/scratch/Software/netMHC-4.0/netMHC"
        netMHC_pan_path = "/homes/users/jtrincado/scratch/Software/netMHCpan-4.0/netMHCpan"
        remove_temp_files = True
        coverage_path = "/homes/users/jtrincado/scratch/test_Junckey/test2/coverageBed/"
        output_path = "/homes/users/jtrincado/scratch/test_Junckey/test2"
        name_user = "jtrincado"
        # ONLY FOR MARVIN
        python2 = "Python/2.7.14-foss-2017b"
        # ONLY FOR HYDRA
        #python2 = "Python/2.7.11"

        #13. Run netMHC-4.0_part2
        logger.info("Part13...")
        run_netMHC_classI_slurm_part2(output_path + "/IR_ORF_filtered_peptide_change.tab", HLAclass_path, HLAtypes_path,
                                      output_path + "/IR_fasta_files",output_path + "/IR_NetMHC-4.0_files", output_path + "/IR_NetMHC-4.0_neoantigens_type_3.tab",
                                      output_path + "/IR_NetMHC-4.0_neoantigens_type_3_all.tab", output_path + "/IR_NetMHC-4.0_neoantigens_type_2.tab",
                                      output_path + "/IR_NetMHC-4.0_neoantigens_type_2_all.tab", output_path + "/IR_NetMHC-4.0_junctions_ORF_neoantigens.tab",
                                      netMHC_path)

        #14. Run netMHCpan-4.0_part2
        logger.info("Part14...")
        run_netMHCpan_classI_slurm_part2(output_path + "/IR_ORF_filtered_peptide_change.tab", HLAclass_path, HLAtypes_pan_path,
                                      output_path + "/IR_fasta_files",output_path + "/IR_NetMHCpan-4.0_files", output_path + "/IR_NetMHCpan-4.0_neoantigens_type_3.tab",
                                      output_path + "/IR_NetMHCpan-4.0_neoantigens_type_3_all.tab", output_path + "/IR_NetMHCpan-4.0_neoantigens_type_2.tab",
                                      output_path + "/IR_NetMHCpan-4.0_neoantigens_type_2_all.tab", output_path + "/IR_NetMHCpan-4.0_junctions_ORF_neoantigens.tab",
                                      netMHC_pan_path)

        logger.info("Done.")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()