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
import cmapPy.math.fast_corr as fast_corr

import plotly.express as pltexpr

logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--verbose", "-v", help="Whether to print a bunch of output.", action="store_true", default=False)
    parser.add_argument("--query_expr_filepath", help="query expression file: sample_id and depmap_id", type=str, required = True)
    parser.add_argument("--ref_expr_filepath", help="reference expression file: in-house RNA-seq dataset", type=str, required = True)
    parser.add_argument("--experiment_id", help = "specific id for experiment", required = True)
    parser.add_argument("--output_subdir", help="subdirectory for output", type=str, default="./cell_line_auth/")
    parser.add_argument("--row_metadata_for_matching", help="the column from query_row_metadata to use to match rows of query with rows of reference", type=str, default = "gene_symbol")
        
    return parser


# reads in the query_expr_file and loads query_expr data
def read_RNA_seq_gctx(query_expr_file):
    # Read in-house RNA seq dataset 
    query_expr_gctoo = parse.parse(query_expr_file)
    query_expr = query_expr_gctoo.data_df
    logger.debug("query_expr.shape \n{}\n".format(query_expr.shape))
    query_expr.head()
    return query_expr_gctoo


# loads the reference expression file into gctoo
def load_ref_expr_data(ref_expr_file):
    ref_expr_gctoo = parse.parse(ref_expr_file)
    ref_expr = ref_expr_gctoo.data_df
    logger.debug("query_expr.shape \n{}\n".format(ref_expr.shape))
    ref_expr.head()
    return ref_expr_gctoo

# adds users choice to query row metadata (Default: gene_index), 
def query_set_index(query_expr_gctoo, row_metadata_for_matching):
    query_expr_df = query_expr_gctoo.data_df.join(
        query_expr_gctoo.row_metadata_df[[row_metadata_for_matching]]
    ).set_index(row_metadata_for_matching)
    return query_expr_df

# returns a subset of the query_expr and subset of the ref_expr that only has the intersection of both indices / shortens size of index (genes)
def build_matched_datasets(query_expr, ref_expr):
    A = set(query_expr.index)
    logger.debug("\nlen(query_expr set):\n{}\n".format(len(A)))
    B = set(ref_expr.index)
    logger.debug("\nlen(ref_expr set):\n{}\n".format(len(B)))

    gene_intersection_set = A.intersection(B)
    logger.debug("\nlen(gene_intersection_set): \n{}\n".format(len(gene_intersection_set)))
    gene_intersection_list = sorted(gene_intersection_set)

    logger.debug("\ngene_intersection_list[:3] : \n{}\n".format(gene_intersection_list[:3]))

    # rearranges the query expression data and reference expression data
    # so that the genes (the index of both dataframes) is in the same order
    subset_query_expr = query_expr.loc[gene_intersection_list]
    subset_ref_expr = ref_expr.loc[gene_intersection_list]

    return subset_query_expr, subset_ref_expr

# fast corr calculation is run
def run_correlation_calculation(subset_query_expr, subset_ref_expr):
    corr_arr = fast_corr.fast_corr(subset_query_expr.to_numpy(), subset_ref_expr.to_numpy())
    logger.debug("\ncorr_arr:\n{}".format(corr_arr))
    corr_df = pd.DataFrame(corr_arr, columns=subset_ref_expr.columns, index = subset_query_expr.columns)
    logger.debug("\ncorr_df\n{}".format(corr_df))
    corr_df.index.name = subset_query_expr.columns.name
    #row sampleid
    #column depmapID
    return corr_df

# save the correlation df as csv/txt file
def save_corr_df(corr_df, exp_id, subdir):
    corr_df_shape = corr_df.shape
    output_filename = "{experiment_id}_cell_line_authentication_corr_r{nrows}x{ncols}.txt".format(
        experiment_id=exp_id, nrows=corr_df_shape[0], ncols=corr_df_shape[1]
    )
    logger.debug("output_filename : {}".format(output_filename))

    output_filepath = os.path.join(subdir, output_filename)
    logger.debug("output_filepath : {}".format(output_filepath))
    
    corr_df.to_csv(output_filepath, sep="\t")
    return output_filepath



def main(args):
    #query expression
    query_expr_gctoo = read_RNA_seq_gctx(args.query_expr_filepath)

    ref_expr_gctoo = load_ref_expr_data(args.ref_expr_filepath)

    query_expr_df = query_set_index(query_expr_gctoo, args.row_metadata_for_matching)
    # query_expr_df = query_expr_gctoo.data_df.join(
    #     query_expr_gctoo.row_metadata_df[[args.row_metadata_for_matching]]
    # ).set_index(args.row_metadata_for_matching)

    subset_query_expr, subset_ref_expr = build_matched_datasets(query_expr_df, ref_expr_gctoo.data_df)

    #run fast_corr calculation
    corr_df = run_correlation_calculation(subset_query_expr, subset_ref_expr) 
 
    # cell_line_auth directory , check if exists
    
    if(os.path.isdir(args.output_subdir)==False):
        os.mkdir(args.output_subdir)

    # save corr df
    output_filepath = save_corr_df(corr_df, args.experiment_id, args.output_subdir)


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)
    logger.debug("args:  {}".format(args))

    main(args)