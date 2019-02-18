"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu
IR_ePydoor.py: get significat exonizations
"""

import os

from lib.IR.extract_significant_IR import *
from lib.IR.IR_associate_gene_ids import *
from lib.IR.IR_kma_associate_gene_ids import *
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

        logger.info("Starting execution IR_ePydoor_part1")

        introns_path = "/projects_rg/SCLC_cohorts/Hugo/Kallisto/iso_tpm_introns.txt"
        bam_path = "/projects_rg/SCLC_cohorts/Hugo/STAR"
        TPM_threshold = 1
        tumor_specific = False
        introns_Normal_path = "/homes/users/jtrincado/scratch/test_Junckey/iso_tpm_introns_Rudin_Normal.txt"
        introns_GTEX_path = "/homes/users/jtrincado/scratch/test_Junckey/chess2.0_assembly_hg19_CrossMap.events_RI_strict.ioe"
        gtf_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.gtf"
        gtf_protein_coding_path = "/projects_rg/SCLC_cohorts/annotation/Homo_sapiens.GRCh37.75.formatted.only_protein_coding.gtf"
        output_path = "/users/genomics/juanluis/SCLC_cohorts/Hugo/epydoor/IR"

        # 0. Format the intron file
        logger.info("Part0...")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        command0 = "module load R; Rscript " + dir_path + "/lib/IR/format_intron_file.R " + introns_path + " " + output_path + "/IR_formatted.tab"
        os.system(command0)

        # 1. Get the IR expressed
        logger.info("Part1...")
        extract_significant_IR(output_path + "/IR_formatted.tab", TPM_threshold, output_path + "/IR_expressed.tab")

        # 2. Obtain the gene ids for the introns.
        logger.info("Part2...")
        # Separate between introns from kma (U2) and U12
        command1="head -n1 "+output_path + "/IR_expressed.tab > "+output_path + "/IR_kma_expressed.tab; grep kma_introns "\
                 +output_path + "/IR_expressed.tab >> "+output_path + "/IR_kma_expressed.tab"
        os.system(command1)
        command2 = "grep -v kma_introns "+output_path + "/IR_expressed.tab > "+output_path + "/IR_no_kma_expressed.tab"
        os.system(command2)
        IR_associate_gene_ids(output_path + "/IR_no_kma_expressed.tab", gtf_path, output_path + "/IR_no_kma_expressed_genes.tab")
        IR_kma_associate_gene_ids(output_path + "/IR_kma_expressed.tab", gtf_path, output_path + "/IR_kma_expressed_genes.tab")
        command3 = "cat "+output_path + "/IR_kma_expressed_genes.tab > "+output_path + "/IR_expressed_genes.tab; tail -n+2 "\
                   +output_path + "/IR_no_kma_expressed_genes.tab >> "+output_path + "/IR_expressed_genes.tab"
        os.system(command3)

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

        # output_path_filtered2 = output_path + "/IR_expressed_genes.tab"
        # 4. Generate random positions for each intron
        logger.info("Part4...")
        generate_random_intronic_positions(output_path_filtered2, gtf_protein_coding_path, 100, output_path + "/random_introns.gtf",
                                           output_path + "/random_introns.bed")

        # 5. Run coverageBed on the samples in the cluster
        logger.info("Part5...")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        #Just for HYDRA
        # command3 = "for sample in $(ls " + bam_path + "/*/*.sorted.bam | cut -d\"/\" -f7 | cut -d\"_\" -f1 | cut -d\".\" -f1 | sort | uniq );do " \
        #                                               "echo \"Processing file $sample: \"$(date); sbatch -J $(echo $sample)_coverageBed " + dir_path + "/coverageBed.sh " + bam_path + "/$(echo $sample)/*.sorted.bam " \
        #                                                    " " + output_path + "/random_introns.bed " + output_path + "/$(echo $sample).coverage_sorted;done"
        #Just for MARVIN
        # command3 = "for sample in $(ls " + bam_path + "/*/*.sorted.bam | cut -d\"/\" -f8 | cut -d\"_\" -f1 | cut -d\".\" -f1 | sort | uniq );do " \
        #                                               "echo \"Processing file $sample: \"$(date); sbatch -J $(echo $sample)_coverageBed " + dir_path + "/coverageBed.sh " + bam_path + "/$(echo $sample)/*.sorted.bam " \
        #                                                                                                                                                                                " " + output_path + "/random_introns.bed " + \
        #            output_path + "/$(echo $sample).coverage_sorted;done"

        # 5.1. Add an extra chrY line ot the end of the -bed file (if don't, coverageBed will crash)
        with open(output_path + "/random_introns.bed", "a") as file:
            file.write("chrY	0	0   Exonization_0_Random_0	+	0")

        # 5.2. Run a job per sample
        command3="for sample in $(ls "+bam_path+"/*/*.bam);do " \
                "sample_id=$(echo $sample | awk -F '/' '{print $(NF-1)}');" \
                "echo \"Processing file $sample: \"$(date); sbatch -J $(echo $sample)_coverageBed "+dir_path+"/coverageBed.sh $(echo $sample) " \
                 + output_path + "/random_introns.bed "+output_path+"/$(echo $sample_id).coverage_sorted;done"
        os.system(command3)
        logger.info("Wait until all jobs have finished. Then, go on with part2")

        exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()