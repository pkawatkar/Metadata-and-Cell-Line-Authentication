from ast import Return
import logging
from operator import not_
from re import T
from shutil import register_unpack_format
import fhtbioinfpy
from fhtbioinfpy.make_genebody_coverage_graphs.make_genebody_coverage_graphs import save_to_tsv
import fhtbioinfpy.setup_logger as setup_logger
import argparse
import sys

import os
from os.path import exists
import pandas as pd
import numpy

import plotly.express as pltexpr

logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--verbose", "-v", help="Whether to print a bunch of output.", action="store_true", default=False)
    parser.add_argument("--metadata_subdir", help="subdirectory for metadata", type=str, default="metadata")
    parser.add_argument("--input_metadata_file_path", "--original_metadata_path", help = "the entire metadata path inputted by user", type = str, required = True)
    parser.add_argument("--experiment_id", help = "specific id for experiment", required = True)
    parser.add_argument("--samples_to_remove_from_metadata", "--samples_to_remove", help = "which samples user desires to remove from metadata, if any.", default=[], nargs="+")
    parser.add_argument("--metadata_columns_to_build_groups", "--selected_cols", help="metadata columns selected to build groups.", default = [], nargs="+")
    return parser


def create_full_path(metadata_subdir, metadata_file):
    input_metadata_dir = metadata_subdir + metadata_file
    logger.debug("Metadata file exists: {}".format(exists(input_metadata_dir)))
    return input_metadata_dir

def load_metadata(input_files):
    orig_metadata = pd.read_csv(
        input_files,
        sep="\t",
        index_col = "sample_id"
    )
    
    #orig_metadata.set_index("sample_id", inplace = True)
    logger.debug("orig_metadata.shape: {}".format(orig_metadata.shape))
    return orig_metadata

def QC_sequence_removal(dataframe, rows_to_remove):
    dataframe = dataframe.drop(rows_to_remove, axis = 1)
    return dataframe


def ask_user_about_replicates():
    user_input = input("Would you like to ignore (I) this error and continue or break (B)?")
    if user_input=="B":
        return True
    
    return False

def check_number_group(metadata_df, metadata_columns_to_build_groups, expected_number_replicates, ask_user_about_replicates=ask_user_about_replicates):
    t = metadata_df[metadata_columns_to_build_groups].value_counts()
    logger.debug('\ndf number of value replicates: \n{}'.format(t))
    not_correct_replicates = t.loc[t!=expected_number_replicates]
    logger.debug('\nvalue replicates not equal to {}\n{}'.format(expected_number_replicates, not_correct_replicates))
    if not_correct_replicates.empty:
        logger.info("Number of value replicates is equal to {}\n".format(expected_number_replicates))
    else:
        logger.info("Number of replicates is not: {}".format(expected_number_replicates))
        do_break = ask_user_about_replicates()

        if do_break:
            raise FhtbioinfpyPrepMetadataCheckValueReplicatesForGroups("Number of replicates is not equal to 3")

    return not_correct_replicates



    #add command line for whether or not user wants to ignore
    #check_number_group(metadata_df, metadata_columns_to_build_groups)

def add_string_columns(dataframe):
    #dataframe['pert_dose_str'] = dataframe['pert_dose'].astype(str) + dataframe['pert_dose_unit']
    #dataframe['pert_time_str'] = dataframe['pert_time'].astype(str) + dataframe['pert_time_unit']

    dataframe = dataframe.assign(pert_dose_str = dataframe['pert_dose'].astype(str) + dataframe['pert_dose_unit'])
    dataframe = dataframe.assign(pert_time_str = dataframe['pert_time'].astype(str) + dataframe['pert_time_unit'])
    return dataframe

def verify(dataframe, group_def_col):
    #locs = all(item in list(dataframe.columns) for item in group_def_col)
    group_def_set = set(group_def_col)
    all_col_set = set(dataframe.columns)
    notFoundSet = group_def_set - all_col_set # tells us which columns were not found in the dataset
    logger.debug('\nset not found:\n{}'.format(notFoundSet))
    #logger.debug('\nCheck that the metadata columns selected to build groups are columns within the metadata file :\n{}'.format(notFoundSet))
    found = False
    if len(notFoundSet) == 0:
        found = True
    elif len(notFoundSet) > 0:
        found = False
        # ___ not found
        print('Some columns selected ({}) were not found in the following set: {}'.format(group_def_set, all_col_set))
        print('\nColumns not found in metadata file:\n{}'.format(notFoundSet))
        
        raise FhtbioinfpyPrepMetadataVerifyColumnsForGroups("Error!")
    #if found == False:
    #    raise FhtbioinfpyPrepMetadataVerifyColumnsForGroups("Error!")

    
    return found

    '''before, verify() returned nothing but I think it's easier if a boolean/locs is returned'''


# takes in the series that you want to add into the original dataframe and cleans it
def clean_selected_series(group_name_series):
    #check this over --> [^A-Za-z0-9]
    #for i in df.columns:
    #print("(")
    group_name_series.astype(str)
    #logger.debug('{}'.format(group_name_series))
    #updated_df = df.copy
    #logger.debug('\nupdated_df:\n{}'.format(updated_df))
    new_series = group_name_series.str.replace("[-)(,]", "")
    new_series = new_series.str.replace(" ", "_")
    logger.debug('\nupdated_series:\n{}'.format(new_series))
    #df = df.str.replace("[-)(,]", "")
    #df = df.str.replace(" ", "_")
    return new_series

    #string that contains

def concatenateGroupDefColumns(cleaned_data, selected_cols):
    selected_cols_df = cleaned_data[selected_cols]
    groups_df = selected_cols_df
    groups_df = groups_df.drop(selected_cols, axis=1) #removal helps because it keeps the number of rows
    logger.debug('\nafter removal:\n{}'.format(groups_df))
    #groups_df['group']= list(pd.Series(groups_df.fillna('').values.tolist()).str.join('_'))
    groups_df['group'] = selected_cols_df.fillna('').apply(lambda row: "_".join(row), axis = 1)
    return groups_df

def add_groups_to_df(original_df, groups_df):
    original_df['group'] = groups_df['group']
    return original_df

def add_experiment_id(df, id):
    df['experiment_id'] = id
    return df


def create_R_Groups(df, col_name):
    # logger.debug('col_name:\n {}'.format(col_name))
    # logger.debug('df:\n {}'.format(df))
    df["R_group"] = df[col_name]
    # logger.debug('df with R_group:\n {}'.format(df))
    # logger.debug('df[col_name]:\n {}'.format(df[col_name]))

    locs = df[col_name].str[0].str.isdigit() # list of whether or not the first char of each item in the column name is a number 
    # logger.debug('locs:\n {}\n'.format(locs))
    
    t = df[col_name][locs]
    #logger.debug('t: \n{}\n'.format(t))
    u = 'x_' + t
    #logger.debug('u:\n {}\n'.format(u))

    df["R_group"][locs] = u
    #logger.debug('df:\n {}\n'.format(df))


    # new_df = df.copy # create copy of df
    # new_df.loc[locs, col_name] = "x_" + new_df.loc[locs, col_name]
    # df["R_group"] = new_df.loc[locs, col_name]
    return df

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
    
    #input_files = 
    #creating metadata directory
    input_files = create_full_path(args.metadata_subdir, args.metadata_file)
    print("input_metadata_dir:  {}".format(input_files))

    #convert files to data frame
    inp_df = load_metadata(input_files)

    #removal of metadata rows not sequenced bc of QC issues
    reduced_df = QC_sequence_removal(inp_df, samples_to_remove_from_metadata)

    #add pert_dose_str and pert_str columns
    metadata_df = add_string_columns(reduced_df)

    #check that "group definition columns" are present in the metadata
    found = verify(metadata_df, metadata_columns_to_build_groups)

    #if found == False:
    #    break

    #remove dashes from columns
    cleaned_data = clean_selected_series(metadata_columns_to_build_groups)

    #concatenate "group definitions" columns to form group names
    df_with_selected_columns = concatenateGroupDefColumns(cleaned_data,  metadata_columns_to_build_groups)

    #add groups to metadata dataframe
    metadata_df = add_groups_to_df(metadata_df, df_with_selected_columns)

    #add experiment id
    metadata_df = add_experiment_id(metadata_df, experiment_id)

    #create R_groups by modifying name of groups to start with an "x_" only if it start with a number
    rGroups_df = create_R_Groups(metadata_df, 'group')

    #save as csv/txt file
    name = "_metadata_r"
        #groups_from_selected_columns_txtfile
    save_to_csv(experiment_id, metadata_df, metadata_subdir, name)
        #R_Groups_txtFile
    save_to_csv(experiment_id, rGroups_df, metadata_subdir, "R_Groups_" + name)

class FhtbioinfpyPrepMetadataVerifyColumnsForGroups(Exception):
    pass

class FhtbioinfpyPrepMetadataCheckValueReplicatesForGroups(Exception):
    pass

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])

    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)