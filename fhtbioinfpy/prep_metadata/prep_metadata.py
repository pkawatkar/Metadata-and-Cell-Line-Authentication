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

import plotly.express as pltexpr

logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--verbose", "-v", help="Whether to print a bunch of output.", action="store_true", default=False)
    parser.add_argument("--metadata_subdir", help="subdirectory for metadata", type=str, default="metadata")
    parser.add_argument("--input_metadata_file_path", help = "the entire metadata path inputted by user", type = str, required = True)
    parser.add_argument("--experiment_id", help = "specific id for experiment", required = True)
    parser.add_argument("--samples_to_remove_from_metadata", help = "which samples user desires to remove from metadata, if any.", default=[], nargs="+")
    parser.add_argument("--metadata_columns_to_build_groups", help="metadata columns selected to build groups.", default = [], nargs="+")
    return parser


def create_full_path(metadata_subdir, metadata_file):
    input_metadata_dir = os.path.join(metadata_subdir, metadata_file)
    logger.debug("Metadata file exists: {}".format(os.path.exists(input_metadata_dir)))
    return input_metadata_dir

def load_metadata(input_file):
    orig_metadata = pd.read_csv(
        input_file,
        sep="\t",
        index_col = "sample_id"
    )
    logger.debug("orig_metadata.shape: {}".format(orig_metadata.shape))
    return orig_metadata

def remove_samples(dataframe, rows_to_remove):
    new_dataframe = dataframe.drop(rows_to_remove, axis = 1)
    return new_dataframe


def ask_user_about_replicates():
    user_input = input("Unexpected number of replicates found for some groups. Would you like to ignore (I) this error and continue or break (B)?")
    return user_input == "B"

def check_replicate_number_per_group(metadata_df, metadata_columns_to_build_groups, expected_number_replicates, ask_user_about_replicates=ask_user_about_replicates):
    group_count_df = metadata_df[metadata_columns_to_build_groups].value_counts()
    logger.debug('\ndf number of value replicates: \n{}'.format(group_count_df))
    unexpected_replicate_count_df = group_count_df.loc[group_count_df!=expected_number_replicates]
    logger.debug('\ngroups with unexpected replicate count - unexpected_replicate_count_df: {}\n{}'.format(expected_number_replicates, unexpected_replicate_count_df))
    if unexpected_replicate_count_df.empty:
        logger.info("All groups have expected number of replicates.\n")
    else:

        do_break = ask_user_about_replicates()
        if do_break:
            msg = """\n\nSome groups had unexpected number of replicates
            expected_number_replicates = {}
            unexpected_replicate_count_df: \n{}\n""".format(expected_number_replicates, unexpected_replicate_count_df)
            logger.exception(msg)
            raise FhtbioinfpyPrepMetadataCheckValueReplicatesForGroups(msg)
    return unexpected_replicate_count_df


def add_string_concated_columns(dataframe):

    dataframe = dataframe.assign(pert_dose_str = dataframe['pert_dose'].astype(str) + dataframe['pert_dose_unit'])
    dataframe = dataframe.assign(pert_time_str = dataframe['pert_time'].astype(str) + dataframe['pert_time_unit'])
    return dataframe

def verify_group_def_columns_in_metadata(dataframe, group_def_col):
    group_def_set = set(group_def_col)
    all_col_set = set(dataframe.columns)
    notFoundSet = group_def_set - all_col_set # tells us which columns were not found in the dataset
    logger.debug('\nset not found:\n{}'.format(notFoundSet))

    if len(notFoundSet) > 0:
        msg = """Some columns selected ({}) were not found in the following set: {} \n
        from the following column: \n{} 
        \nColumns not found in metadata file:\n{}""".format(group_def_set, all_col_set, group_def_col, notFoundSet)
        logger.exception(msg)
        raise FhtbioinfpyPrepMetadataVerifyColumnsForGroups(msg)

    elif len(notFoundSet) == 0:
        return True


# takes in the series that you want to add into the original dataframe and cleans it
def clean_selected_series(group_name_series):

    group_name_series = group_name_series.astype(str)
    new_series = group_name_series.str.replace("[-)(,]", "")
    logger.debug('\nupdate1:\n{}'.format(new_series))
    new_series = new_series.str.replace(" ", "_")
    logger.debug('\nupdate2:\n{}'.format(new_series))
    new_series = new_series.str.replace("[^A-Za-z0-9_]+","")
    logger.debug('\nupdate3:\n{}'.format(new_series))
    #logger.debug('\nupdated_series:\n{}'.format(new_series))
    return new_series


def build_group_names(cleaned_metadata_df, group_def_cols):
    group_def_cols = cleaned_metadata_df[group_def_cols]
    groups_df = group_def_cols
    groups_df = groups_df.drop(group_def_cols, axis=1) #removal helps because it keeps the number of rows
    logger.debug('\nafter removal:\n{}'.format(groups_df))
    groups_df['group'] = group_def_cols.fillna('').apply(lambda row: "_".join(row), axis = 1)
    return groups_df

def add_groups_to_df(original_df, groups_series):
    original_df['group'] = groups_series
    return original_df

def add_experiment_id(df, experiment_id):
    df['experiment_id'] = experiment_id
    return df


def create_R_Groups(original_df, col_name):
    original_df["R_group"] = original_df[col_name]
    logger.debug('df with R_group:\n {}'.format(original_df))
    logger.debug('df[col_name]:\n {}'.format(original_df[col_name]))
    locs = original_df[col_name].str[0].str.isdigit() # list of whether or not the first char of each item in the column name is a number 
    logger.debug('locs:\n {}\n'.format(locs))

    new_df = original_df[col_name][locs]
    logger.debug('new_df: \n{}\n'.format(new_df))
    new_df_withX = 'x_' + new_df
    logger.debug('new_df_withX:\n {}\n'.format(new_df_withX))
    original_df["R_group"][locs] = new_df_withX
    logger.debug('df:\n {}\n'.format(original_df))

    return original_df

def build_output_file_path(id, metadata_subdir, name, df_shape):
    output_filename = id + name + str(df_shape[0]) + 'x' + str(df_shape[1]) + '.txt'
    logger.debug("output_filename : {}".format(output_filename))

    output_filepath = os.path.join(metadata_subdir, output_filename)
    logger.debug("output_filepath : {}".format(output_filepath))

    return output_filepath



def save_to_csv(output_filepath, df):
    df.to_csv(output_filepath, sep="\t")


def main(args):
    experiment_id = "NS-21.0138"
    samples_to_remove_from_metadata = [] # ["G8M54"]
    metadata_subdir = "./metadata/"
    metadata_file = "2021-11-09 next-gen-sequencing-annotations-ATACseq-KA.txt"
    metadata_columns_to_build_groups = ["pert_id"]
 
    #creating metadata directory
    input_files = create_full_path(args.metadata_subdir, args.metadata_file)
    print("input_metadata_dir:  {}".format(input_files))

    #convert files to data frame
    inp_df = load_metadata(input_files)

    #removal of metadata rows not sequenced bc of QC issues
    reduced_df = remove_samples(inp_df, samples_to_remove_from_metadata)

    #add pert_dose_str and pert_str columns
    metadata_df = add_string_concated_columns(reduced_df)

    #check that "group definition columns" are present in the metadata
    found = verify_group_def_columns_in_metadata(metadata_df, metadata_columns_to_build_groups)

    #remove dashes from columns
    cleaned_data = clean_selected_series(metadata_columns_to_build_groups)

    #concatenate "group definitions" columns to form group names
    df_with_selected_columns = build_group_names(cleaned_data,  metadata_columns_to_build_groups)

    #add groups to metadata dataframe
    metadata_df = add_groups_to_df(metadata_df, df_with_selected_columns)

    #add experiment id
    metadata_df = add_experiment_id(metadata_df, experiment_id)

    #create R_groups by modifying name of groups to start with an "x_" only if it start with a number
    rGroups_df = create_R_Groups(metadata_df, 'group')

    #save as csv/txt file
    name = "_metadata_r"
    # build output file path
    groups_from_selected_columns_output_filepath = build_output_file_path(experiment_id, metadata_subdir, name, metadata_df.shape)
    rGroups_output_filepath = build_output_file_path(experiment_id, metadata_subdir, "R_Groups_" + name, rGroups_df.shape)
        #groups_from_selected_columns_txtfile
    save_to_csv(groups_from_selected_columns_output_filepath, metadata_df)
        #R_Groups_txtFile
    save_to_csv(rGroups_output_filepath, rGroups_df)

class FhtbioinfpyPrepMetadataVerifyColumnsForGroups(Exception):
    pass

class FhtbioinfpyPrepMetadataCheckValueReplicatesForGroups(Exception):
    pass

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)
    logger.debug("args:  {}".format(args))

    main(args)