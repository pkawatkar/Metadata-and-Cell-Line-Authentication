""" 
 - code looks up the depmap id using the biocontextid -> adds a column -> in prep metadata(cell-line-lookup
 - later on steps: 
      - after the differential analysis , cell line Authentication
      - only will need to run the correlation calculation 
"""
import ast
import logging
import operator
import re
import shutil
import fhtbioinfpy.setup_logger as setup_logger
import argparse
import sys

import os
import pandas as pd
import numpy
import collections
import cmapPy
import cmapPy.pandasGEXpress.parse as parse

import plotly.express as pltexpr

logger = logging.getLogger(setup_logger.LOGGER_NAME)

def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--verbose", "-v", help="Whether to print a bunch of output.", action="store_true", default=False)
    return parser


def add_cleaned_cl_name(metadata, col_name):
    #type of col_name should be string assertTrue(type(col_name)=="str")
    metadata["stripped_" + col_name] = metadata[col_name]
    logger.debug("\nmetadata\n{}".format(metadata))
    # 1) force all names to be uppercase 
    metadata["stripped_" + col_name] = metadata["stripped_" + col_name].astype(str).str.upper()
    # 2) remove all special characters before comparing
    metadata["stripped_" + col_name] = metadata["stripped_" + col_name].astype(str).str.replace('[^A-Za-z0-9]+','')
    # 3) remove LabGuru_cleaned_annot_cl
    metadata["stripped_" + col_name] = metadata["stripped_" + col_name].astype(str).str.replace('NCI','')
    return metadata

def shorten_sample_info(sample_info_file):
    sample_info = pd.read_csv(sample_info_file)
    # sample_info.set_index(['DepMap_ID'], inplace=True)
    sample_info = sample_info[['DepMap_ID', 'stripped_cell_line_name']]
    return sample_info

#metadata has stripped_bio_context_id,
#sample_info has stripped_cell_line_name
#exception: 1) stripped_bio_context_id is never equal to a stripped_cell_line_name and no corresponding depmapID is found

def verify_match(metadata, sample_info):
    look_up_df = sample_info.set_index("stripped_cell_line_name")
    #change index in sample info to stripped_cell_line_name
    joined_dataframe = metadata.join(look_up_df, on = "stripped_bio_context_id", how = "left")
    logger.debug("\njoined_dataframe:\n\n{}".format(joined_dataframe))
    isNull_df = joined_dataframe["DepMap_ID"].isnull()
    logger.debug("\nisNull_df:\n\n{}".format(isNull_df))
    dmID_noVal = joined_dataframe.loc[isNull_df]
    logger.debug("\nno match values:\n\n{}".format(dmID_noVal))
    if not dmID_noVal.empty:
        msg = """\n\n!!!\nSome stripped_bio_context_id values don't have an equal stripped_cell_line_name.\n No corresponding DepMap_ID found.\n\n{}""".format(dmID_noVal)
        logger.exception(msg)
        raise FhtbioinfpydepmapIDVerifyEquivalenceOfCLMandBCID(msg)

    return joined_dataframe


'''
def verify_match(metadata, sample_info):
    for i in metadata["stripped_bio_context_id"]:
        neverEqual = False
        for j in sample_info["stripped_cell_line_name"]:
            if i == j:
                neverEqual = True
                index_val = sample_info[sample_info["stripped_cell_line_name"]==i].index.values #find the index value of specific stripped_cell_line_name
                depmapID = sample_info.loc[index_val, "DepMap_ID"] # find the corresponding depmapID at that index value
                if depmapID.isnan():
                    msg = """\n\n!!!\nCorresponding DepMap_ID not found\n metadata[stripped_bio_context_id]: \n{}\nsample_info[stripped_cell_line_name]: \n{}\n""".format(i, j)
                    logger.exception(msg)
                    raise FhtbioinfpydepmapIDCorrespondingDepmapIDNotFound(msg)
        if neverEqual == False:
            msg = """\n\n!!!\nStripped_bio_context_id (  {}  ) is never equal to a stripped_cell_line_name\n\n""".format(i)
            logger.exception(msg)
            raise FhtbioinfpydepmapIDVerifyEquivalenceOfCLMandBCID(msg)

    return neverEqual
'''

class FhtbioinfpydepmapIDVerifyEquivalenceOfCLMandBCID(Exception):
    pass
