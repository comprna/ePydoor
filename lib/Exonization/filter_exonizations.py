"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

filter_exonizations: from our list opf exonizations, we are gonna remove those ones that appears in the Rudin
and Intropolis data. We will exclude a junction if we see that there is at least a sample in any data that has a
minimum number of reads
"""

import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import logging, sys, os, re
from collections import Counter

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

        # exonizations_path = sys.argv[1]
        # rudin_path = sys.argv[2]
        # intropolis_path = sys.argv[3]
        # output_path = sys.argv[4]
        # flag_Rudin = sys.argv[5]


        # ## exonizations_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/exonizations_by_sample_coverage.tab"
        # exonizations_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/non_mutated_exonizations.tab"
        # rudin_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering/exonizations_by_sample_Rudin_normal_junctions.tab"
        # intropolis_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/exonizations_Intropolis_by_sample.tab"
        # output_path = "/projects_rg/SCLC_cohorts/George/PSI_Junction_Clustering_v2/non_mutated_exonizations_by_sample_filtered.tab"
        # flag_Rudin = "Rudin"

        if(flag_Rudin!="No_Rudin"):

            #Load the exonizations
            exonizations = pd.read_table(exonizations_path, delimiter="\t")
            #Load the Rudin exons
            rudin = pd.read_table(rudin_path, delimiter="\t")
            #Load the Intropolis exons
            intropolis = pd.read_table(intropolis_path, delimiter="\t")

            #Merge the files and extract the exons that are not in the other files
            merge1 = exonizations.merge(rudin, on='New_exon', how='left', suffixes=('', '_y'))
            merge1_f = merge1[~merge1.Sample_id_y.notnull()]
            #Get also the samples that would be filtered because of Rudin
            merge1_f_null = merge1[merge1.Sample_id_y.notnull()]
            sample_counts = merge1_f_null["Sample_id"].value_counts()
            #Save this in an external file
            path1 = "/".join(output_path.split("/")[:-1])
            sample_counts.to_csv(path1+"/samples_removed_Rudin.tab", sep="\t", index=True, header=True)
            merge2 = merge1_f.merge(intropolis, on='New_exon', how='left', suffixes=('', '_z'))
            merge2_f = merge2[~merge2.Sample_id_z.notnull()]

            #Remove the extra columns
            df = merge2_f.iloc[:, 0:16].copy()
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
            merge1 = exonizations.merge(intropolis, on='New_exon', how='left', suffixes=('', '_y'))
            merge1_f = merge1[~merge1.Sample_id_y.notnull()]

            # Remove the extra columns
            df = merge1_f.iloc[:, 0:12].copy()
            df.to_csv(output_path, sep="\t", index=False, header=True)
            logger.info("Saved " + output_path)
            logger.info("Done. Exiting program.")
            exit(0)


    except Exception as error:
        logger.error('ERROR: ' + repr(error))
        logger.error("Aborting execution")
        sys.exit(1)
