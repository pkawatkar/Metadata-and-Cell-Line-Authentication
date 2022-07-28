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
    # parser.add_argument("--ref_expr_filepath", help="reference expression file: in-house RNA-seq dataset", type=str, required = True)
    parser.add_argument("--experiment_id", help = "specific id for experiment", required = True)
    parser.add_argument("--output_subdir", help="subdirectory for output", type=str, default="./cell_line_auth/")
    parser.add_argument("--corr_filepath", help="correlation_depmap filepath to load the correlation df", type=str, required = True)        
    return parser

#load the correlation df that was saved in cell_line_correlation
def load_correlation_df(corr_filepath):
    corr_df = pd.read_csv(corr_filepath, sep='\t', index_col = 0)
    logger.debug("\ncorr_df\n{}".format(corr_df))    
    # corr_df.set_index("cid", inplace = True)
    logger.debug("\ncorr_df\n{}".format(corr_df))
    return corr_df

# reload query_expr gctoo
def read_RNA_seq_gctx(query_expr_file):
    # Read in-house RNA seq dataset 
    query_expr_gctoo = parse.parse(query_expr_file)
    query_expr = query_expr_gctoo.data_df
    logger.debug("query_expr.shape \n{}\n".format(query_expr.shape))
    query_expr.head()
    return query_expr_gctoo

# return series of the depmap ID with the best correlation
def identify_rank(corr_rank_df):
    cl_match_ID_s = corr_rank_df.idxmin(axis=1) #depmap id with best correlation
    logger.debug("\ncl_match_ID_s\n{}".format(cl_match_ID_s))
    cl_match_ID_s.name = "top_corr_depmap_ID"
    return cl_match_ID_s

# adds column onto df that tells us whether the depmap_id that has the best correlation matches the depmap id in query_expr col metadata
def compare_depmap_id(query_expr_col_meta_df, cl_match_ID_s):
    logger.debug("\n\nquery_expr_col_meta_df\n {} \n df_cl_match_ID_depmap_ID_col\n{}\n".format(query_expr_col_meta_df, cl_match_ID_s))
    df_match = query_expr_col_meta_df.join(cl_match_ID_s)
    logger.debug("\n df_match\n{}".format(df_match))
    df_match["cell_line_match"] = df_match.DepMap_ID == df_match["top_corr_depmap_ID"]
    logger.debug("\n df_match\n{}".format(df_match))
    return df_match
    
# checks for whether any values in the cell_line_match column of df are False, meaning the depmap_id doesn't have a match
def check_for_warnings(compared_depmaps_df):
    #depmapid, bio_context_id, top_corr, cell_linematch
    depmap_no_match_df = compared_depmaps_df[~compared_depmaps_df.cell_line_match]
    logger.debug("\n~compared_depmaps_df.cell_line_match:\n{}".format(~compared_depmaps_df.cell_line_match))
    logger.debug("\ndepmap_no_match_df:\n{}".format(depmap_no_match_df))
    if not depmap_no_match_df.empty:
        msg = """\n!!!\n The following DepMap IDs have no match. \n{}\n""".format(depmap_no_match_df)
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):
            logger.exception(msg)
            print(depmap_no_match_df[["DepMap_ID", "bio_context_id", "top_corr_depmap_ID", "cell_line_match"]])
            raise FhtbioinfpyCellLineAuthenticationNoDepMapIDMatch(msg)
     

# saves the compared_depmap_df as a csv/txt file
def save_compared_depmap_df(compared_depmap_df, exp_id, subdir):
    compared_depmap_df_shape = compared_depmap_df.shape
    output_filename = "{experiment_id}_cell_line_authentication_compared_depmap_r{nrows}x{ncols}.txt".format(
        experiment_id=exp_id, nrows=compared_depmap_df_shape[0], ncols=compared_depmap_df_shape[1]
    )
    logger.debug("output_filename : {}".format(output_filename))
    output_filepath = os.path.join(subdir, output_filename)
    logger.debug("output_filepath : {}".format(output_filepath))
    compared_depmap_df.to_csv(output_filepath, sep="\t")
    return output_filepath

def main(args):
    query_expr_gctoo = read_RNA_seq_gctx(args.query_expr_filepath)
    # load correlation_df
    corr_df = load_correlation_df(args.corr_filepath)
    
    if not os.path.isdir(args.output_subdir):
        os.mkdir(args.output_subdir)

   
    corr_rank_df = corr_df.rank(axis = 1, method = "first", ascending = False)
    cl_match_ID_s = identify_rank(corr_rank_df) #1 of the columns


    #yes/no column - does query_exp depmap match with ref_Exp depmap
    compared_depmap_df = compare_depmap_id(query_expr_gctoo.col_metadata_df, cl_match_ID_s) #merge them and check if they match

    # save df
    save_compared_depmap_df(compared_depmap_df, args.experiment_id, args.output_subdir)

    #looks through output_df, finds where compared)depmaps_df is False , throw exception "doesn't match"
    check_for_warnings(compared_depmap_df)





class FhtbioinfpyCellLineAuthenticationNoDepMapIDMatch(Exception):

    pass


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)
    logger.debug("args:  {}".format(args))

    main(args)

