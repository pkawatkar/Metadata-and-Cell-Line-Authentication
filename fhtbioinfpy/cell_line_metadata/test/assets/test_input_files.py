import os.path
import cmapPy.pandasGEXpress.parse as parse
import cmapPy.pandasGEXpress.subset_gctoo as sg
import cmapPy.pandasGEXpress.write_gct as write_gct

import fhtbioinfpy.cell_line_metadata.cell_line_metadata as clm



def subset_input_files(query_expr_file):
    
    query_expr_df = parse.parse(query_expr_file) #load query expression from filepath
    print(query_expr_df)
    query_expr_subset_df = sg.subset_gctoo(query_expr_df, row_bool=None, col_bool=None, rid=None, cid=None,
                ridx=list(range(1000)), cidx=None, exclude_rid=None, exclude_cid=None) #selects first 1000 rows
    print(query_expr_subset_df)
    write_gct.write(query_expr_subset_df, "./test_example_input_data/query_expr_subset.gct") #write and save
    #write_gct.write(query_expr_subset_df, "./query_expr_subset.gct") #write and save

    return query_expr_subset_df


if __name__ == "__main__":
    query_expr_file = "./example_input_data/NS-22.0019_log2_TPM_r57771x18.gct" #human genes
    subset_input_files(query_expr_file)