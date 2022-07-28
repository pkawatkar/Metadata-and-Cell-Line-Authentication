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
    parser.add_argument("--query_expr_file", help="query expression file: in-house RNA-seq dataset", type=str, required = True)
    parser.add_argument("--sample_info_file", help = "sample info containing DepMap_ID and the cell line name", type = str, required = True)
    
    # # In-house RNA-seq dataset
    # query_expr_file = "./input_data/2022-02-10_heatmap_log2_TPM_r113153x1115.gctx"
    # # Read the sample_info containing the DepMap_ID and the cell line name
    # sample_info_file = "./input_data/sample_info.csv"
        
    return parser



def read_RNA_seq_gctx(query_expr_file):
    # Read in-house RNA seq dataset 
    query_expr_gctx = parse.parse(query_expr_file)
    query_expr = query_expr_gctx.data_df

    print(query_expr.shape)
    query_expr.head()
    return query_expr_gctx, query_expr


def remove_mouse_genes(query_expr):
    new_dfq = query_expr[query_expr.index.str.contains("ENSMUSG") == False]
    logger.debug("query without mouse genes: \n {}".format(new_dfq))
    return new_dfq

def drop_nas_query_df(dataframe):
    '''remove genes that have some nas (not all)?'''
    # Check for NAs in index
    NAs = dataframe.index.isna()
    indexNAs = collections.Counter(list(NAs))
    logger.debug("NAs in the indexes: \n {}".format(indexNAs))
    # Check for NAs in columns
    NAs = dataframe.columns.isna()
    columnNAs = collections.Counter(list(NAs))
    logger.debug("NAs in the columns: \n {}".format(columnNAs))
    # Check for NaN in DataFrame
    NAs = dataframe.isnull().all()
    dataframeNAs = collections.Counter(list(NAs))
    logger.debug("NAs in the dataframe: \n {}".format(dataframeNAs))

    # View columns with all NAs (Mouse samples) if any
    if dataframeNAs[True] == 0:
        print("There is no columns with only NAs in the Dataframe.")
        NA_samples = []
    else:
        print('There is', dataframeNAs[True] , 'columns with only NAs in this Dataframe.')
        NA_samples = list(dataframe[dataframe.columns[dataframe.isnull().all()]].columns)
        print('Names of columns/samples with all NAs in DataFrame: \n', NA_samples)
        # Drop columns with all NAs in query_expr
        dataframe = dataframe.dropna(axis='columns')
        assert(collections.Counter(list(dataframe.isnull().all()))[True]==0)
    return dataframe, NA_samples


def read_metadata_col_gctx(query_expr_gctx, NA_samples):
    # Read the metadata from 10_heatmap_log2_TPM_r113153x1115.gctx to compare the cell line sample 
    # and the cell line with the highest correlation
    metadata = pd.DataFrame(query_expr_gctx.col_metadata_df['bio_context_id'])
    metadata.index.name = "sample"
    metadata.drop(NA_samples, inplace= True) # Drop the samples with all NAs samples (Mouse samples)
    metadata.reset_index(inplace=True)
    return metadata


def drop_nas_metadata(metadata):
    # Check for NAs in index
    NAs = metadata.index.isna()
    indexNAs = collections.Counter(list(NAs))
    logger.debug("NAs in the indexes: \n {}".format(indexNAs))

    # Check for NAs in columns
    NAs = metadata.columns.isna()
    columnNAs = collections.Counter(list(NAs))
    logger.debug("NAs in the columns: \n {}".format(columnNAs))

    # Check for NaN in DataFrame
    NAs = metadata.isnull().all()
    dataframeNAs = collections.Counter(list(NAs))
    logger.debug("NAs in the dataframe: \n {}".format(dataframeNAs))

    # View samples with no annotation if any
    sample_wth_no_annot = [sample for sample, row in metadata.iterrows() if row.isnull().any()]
    if sample_wth_no_annot==[]:
        print('All samples have an annotation.')
    else:
        print('Samples with no annotation: \n', sample_wth_no_annot)
        # Replace samples with no annotation by 'NoAnnotation' in metadata if NAs
        metadata.fillna('NOANNOTATION', inplace=True)

    # Sanity Check
    print('New # of NAs in DataFrame: ', metadata.isnull().all())
    print('Dim of the DataFrame: ', metadata.shape)
    return metadata

#following method used to clean the metadata and the sample info
'''
def add_cleaned_cl_name(metadata, col_name):
    #type of col_name should be string assertTrue(type(col_name)=="str")
    metadata["cleaned_" + col_name] = metadata[col_name]
    logger.debug("\nmetadata\n{}".format(metadata))
    # 1) force all names to be uppercase 
    metadata["cleaned_" + col_name] = metadata["cleaned_" + col_name].astype(str).str.upper()
    # 2) remove all special characters before comparing
    metadata["cleaned_" + col_name] = metadata["cleaned_" + col_name].astype(str).str.replace('[^A-Za-z0-9]+','')
    # 3) remove LabGuru_cleaned_annot_cl
    metadata["cleaned_" + col_name] = metadata["cleaned_" + col_name].astype(str).str.replace('NCI','')
    return metadata
'''
'''
def shorten_sample_info(sample_info_file):
    sample_info = pd.read_csv(sample_info_file)
    # sample_info.set_index(['DepMap_ID'], inplace=True)
    sample_info = sample_info[['DepMap_ID', 'stripped_cell_line_name']]
    return sample_info
'''
def drop_nas_sample_info(sample_info):
    # Check for NAs in index
    logger.debug("\nNAs in the DepMap_ID indexes: {}".format(collections.Counter(sample_info.index.isna())))
    # Check for NAs in columns
    logger.debug("\nNAs in the columns: {}".format(collections.Counter(list(sample_info.columns.isna()))))
    # Check for NaN in DataFrame
    logger.debug("\nNAs in DataFrame: {}".format(sample_info.isnull().all()))
    # Subset rows with NAs
    logger.debug("\nRow(s) with NA: \n{}".format(sample_info[sample_info.isnull().any(axis=1)]))
    # Drop NAs
    sample_info.dropna(subset = ["stripped_cell_line_name"], inplace=True)
    # Sanity Check
    logger.debug("\nNew # of NAs in DataFrame: {}".format(sample_info.isnull().all()))
    logger.debug("\nDim of the DataFrame: \n{}".format(sample_info.shape[0]))

    return sample_info

def merge_metadata_w_sample_info(metadata, sample_info):

    A = set(metadata["cleaned_bio_context_id"])

    B= set(sample_info["cleaned_stripped_cell_line_name"])
    isEqual =  A.issubset(B)
    # set(metadata["cleaned_bio_context_id"]) == sample_info["cleaned_stripped_cell_line_name"]
    logger.debug("isEqual: {}".format(isEqual))
    # logger.debug("isEqual.all: {}".format(isEqual.all))
    if not isEqual:
    #if metadata["cleaned_bio_context_id"].ne(sample_info["cleaned_stripped_cell_line_name"]):
        msg = """Unable to merge data.\nThe cleaned_bio_context_id column of the metadata and the cleaned_stripped_cell_line_name column of the sample_info are not equal.\n
        cleaned_bio_context_id column:  {}\n
        cleaned_stripped_cell_line_name column:  {}\n
        """.format(metadata["cleaned_bio_context_id"], sample_info["cleaned_stripped_cell_line_name"])
        logger.exception(msg)
        raise FhtbioinfpyCellLineMetadataCheckEquivalenceOfBCIDandStrippedCLName(msg)
    sample_df = metadata.merge(sample_info, 
                                how='left', 
                                left_on='cleaned_bio_context_id', 
                                right_on='cleaned_stripped_cell_line_name').drop('cleaned_stripped_cell_line_name', 1)
    logger.debug("\nsample_df.shape{}".format(sample_df.shape))
    logger.debug("\nsample_df.tail(): \n {}".format(sample_df.tail()))
    return sample_df



def build_output_file_name(df_shape):
    output_filename = "cell_line_metadata{nrows}x{ncols}.txt".format(
        nrows=df_shape[0], ncols=df_shape[1]
    )
    logger.debug("output_filename : {}".format(output_filename))

    return output_filename

def main(args):
    #query expression
    query_expr_gctx, query_expr = read_RNA_seq_gctx(args.query_expr_file)
    dataframe_query = remove_mouse_genes(query_expr) #not necessary
    dataframe_query, NA_samples = drop_nas_query_df(dataframe_query) #take look into
    
    output_filename = build_output_file_name(dataframe_query.shape)
    subdir = "assets/DepMapIDmetadata/"
    output_filepath = os.path.join(subdir, output_filename)

    #SAVE DATAFRAME QUERY HERE
    dataframe_query.to_csv(output_filepath, sep="\t")

    #metadata
    metadata = read_metadata_col_gctx(query_expr_gctx, NA_samples)
    metadata = drop_nas_metadata(metadata) #take look into
   # metadata = add_cleaned_cl_name(metadata, "bio_context_id")

    #sample info
    #sample_info = shorten_sample_info(args.sample_info_file)
    sample_info = drop_nas_sample_info(args.sample_info_file) #take look into
    #sample_info = add_cleaned_cl_name(sample_info, "stripped_cell_line_name")
    
    #merging metadata and sample info
    merged_data = merge_metadata_w_sample_info(metadata, sample_info)

    output_filename = build_output_file_name(merged_data.shape)
    subdir = "assets/DepMapIDmetadata/"
    output_filepath = os.path.join(subdir, output_filename)

    #SAVE MERGED DATA HERE
    merged_data.to_csv(output_filepath, sep="\t")


class FhtbioinfpyCellLineMetadataCheckEquivalenceOfBCIDandStrippedCLName(Exception):
    pass

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)
    logger.debug("args:  {}".format(args))

    main(args)