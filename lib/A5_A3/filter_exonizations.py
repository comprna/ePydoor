"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

filter_exonizations: from our list opf exonizations, we are gonna remove those ones that appears in the other datsets
like Rudin or Intropolis
"""

import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import logging, sys, os, re

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


def filter_exonizations(exonizations_path, rudin_path, intropolis_path, output_path, flag_Rudin):

    try:
        logger.info("Starting execution")
        #
        ## exonizations_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/new_A5_A3_junctions_by_sample.tab"
        # exonizations_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/non_mutated_new_A5_A3_junctions.tab"
        # rudin_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/new_A5_A3_junctions_by_sample_Rudin_normal.tab"
        # intropolis_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/new_A5_A3_junctions_by_sample_Intropolis.tab"
        # output_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/non_mutated_new_A5_A3_junctions_filtered.tab"
        # flag_Rudin = "Rudin"

        if(flag_Rudin!="No_Rudin"):

            #Load the eoxnizations
            exonizations = pd.read_table(exonizations_path, delimiter="\t", )
            #Load the Rudin exons
            rudin = pd.read_table(rudin_path, delimiter="\t")
            #Load the Intropolis exons
            intropolis = pd.read_table(intropolis_path, delimiter="\t")

            #Merge the files and extract the exons that are not in the other files
            merge1 = exonizations.merge(rudin, on=['Canonical_Junction_id','Alt_Junction_id'], how='left', suffixes=('', '_y'))
            merge1_f = merge1[~merge1.Gene_y.notnull()]
            merge2 = merge1_f.merge(intropolis, on=['Canonical_Junction_id','Alt_Junction_id'], how='left', suffixes=('', '_z'))
            merge2_f = merge2[~merge2.Gene_z.notnull()]

            #Remove the extra columns
            df = merge2_f.iloc[:, 0:19].copy()
            # df = merge2_f.copy()
            df.to_csv(output_path, sep="\t", index=False, header=True)
            logger.info("Saved "+output_path)
            logger.info("Done. Exiting program.")
            exit(0)

        else:

            # Load the exonizations
            exonizations = pd.read_table(exonizations_path, delimiter="\t")
            # Load the Intropolis exons
            intropolis = pd.read_table(intropolis_path, delimiter="\t")

            # Merge the files and extract the exons that are not in the other files
            merge1 = exonizations.merge(intropolis, on=['Canonical_Junction_id','Alt_Junction_id'], how='left', suffixes=('', '_y'))
            merge1_f = merge1[~merge1.Sample_id_y.notnull()]

            # Remove the extra columns
            df = merge1_f.iloc[:, 0:16].copy()
            # df = merge1_f.copy()
            df.to_csv(output_path, sep="\t", index=False, header=True)
            logger.info("Saved " + output_path)
            logger.info("Done. Exiting program.")
            exit(0)

    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
