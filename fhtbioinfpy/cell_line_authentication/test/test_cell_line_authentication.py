from cgi import test
from re import A
import unittest
import logging
import fhtbioinfpy.setup_logger as setup_logger
import fhtbioinfpy.cell_line_authentication.cell_line_authentication as cla
import pandas as pd
import os
import tempfile
import numpy as np
import cmapPy
import cmapPy.pandasGEXpress.parse as parse

logger = logging.getLogger(setup_logger.LOGGER_NAME)

class TestCellLineAuthentication(unittest.TestCase):

    def setUp(self):
        logger.debug("setUp")   
        
        self.query_expr_filepath = "./test_cell_line_auth/test_cell_line_auth_query_data_r110x86.gctx"
        self.ref_expr_filepath   = "./test_cell_line_auth/test_cell_line_auth_ref_data_r108x1406.gctx"

        self.df = pd.DataFrame({"c":[-1.0, 1.0], "d":[1.0, -1.0], "e":[1.0, -1.0]})
        self.df = self.df.set_index(pd.Series(['a', 'b']))

        self.query_expr_gctoo = parse.parse(self.query_expr_filepath)
        logger.debug('\nself.df.shape\n{}'.format(self.df.shape))
        logger.debug('\nself.df.head\n{}'.format(self.df.head()))
        
        self.exp_id = "test_experiment_id"
        self.subdir = "./assets/cell_line_auth/"
        self.corr_filepath = "./test_cell_line_auth/test_exp_id_cell_line_authentication_corr_r86x1406.txt"
        

    def tearDown(self):
        logger.debug("tearDown")


    def test_main_functional(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            logger.debug("test_main tmpdirname:  {}".format(tmpdirname))
            
            arg_list = ["--query_expr_filepath", self.query_expr_filepath,
                        "--experiment_id", "test_experiment_id",
                        "--output_subdir", tmpdirname,
                        "--corr_filepath", self.corr_filepath

            ]

            args = cla.build_parser().parse_args(arg_list)
            logger.debug(args)
            with self.assertRaises(cla.FhtbioinfpyCellLineAuthenticationNoDepMapIDMatch) as context:
                cla.main(args)
            expected_output_file = os.path.join(tmpdirname, "test_experiment_id_cell_line_authentication_compared_depmap_r86x35.txt") # from build_output_file_name
            logger.debug("expected_output_file: {}".format(expected_output_file))
            self.assertTrue(os.path.isfile(expected_output_file))

            loaded_csv_file = pd.read_csv(expected_output_file, sep = "\t")
            self.assertFalse(loaded_csv_file.empty)

    def test_load_correlation_df(self):
        logger.debug("\ntest_read_RNA_seq_gctx\n")

        self.assertTrue(os.path.exists(self.corr_filepath))
        # check the two dataframes are equal
        corr_df = cla.load_correlation_df(self.corr_filepath)
        
        self.assertFalse(corr_df.empty)
        logger.debug("\ncorr_df.shape\n{}".format(corr_df.shape))
        self.assertEqual(corr_df.shape, (86, 1406)) 
        logger.debug("\nquery_expr_df.iloc[1,2]:\n{}".format(corr_df.iloc[1,2])) 
        self.assertEqual(corr_df.iloc[1,2], 0.7613657886411502)
    
    def test_read_RNA_seq_gctx(self):
        logger.debug("\ntest_read_RNA_seq_gctx\n")

        self.assertTrue(os.path.exists(self.query_expr_filepath))
        # check the two dataframes are equal
        query_expr_gctoo = cla.read_RNA_seq_gctx(self.query_expr_filepath)
        query_expr_df = query_expr_gctoo.data_df

        logger.debug("\nquery_expr_data.shape: {}".format(query_expr_df.shape))   
        self.assertFalse(query_expr_df.empty)
        self.assertEqual(query_expr_df.shape, (110, 86)) 
        logger.debug("\nquery_expr_df.iloc[1,2]:\n{}".format(query_expr_df.iloc[1,2])) 
        self.assertEqual(query_expr_df.iloc[1,2], 3.307300090789795) 


# corr_rank_df
    def test_identify_rank(self):
        logger.debug("\ntest_identify_rank\n")

        # edge case 1: all ranks are different, no repeats
        rank_df = pd.DataFrame({"depmap_id1":[1, 2, 3], "depmap_id2":[3, 1, 2], "depmap_id3":[2, 3, 1]})
        rank_df = rank_df.set_index(pd.Series(['a', 'b', 'c']))
        logger.debug("\nrank_df\n{}".format(rank_df))
        cl_match_ID_s = cla.identify_rank(rank_df)
        logger.debug("\ncl_match_ID_s\n{}".format(cl_match_ID_s))
        expected_series = pd.Series(["depmap_id1", "depmap_id2", "depmap_id3"], index = ['a', 'b', 'c'])
        logger.debug("\nexpected_series\n{}".format(expected_series))

        isEqual = expected_series == cl_match_ID_s
        logger.debug("\nisEqual:\n{}".format(isEqual))
        self.assertTrue(isEqual.all())


    def test_compare_depmap_id(self):
        logger.debug("\ntest_compare_depmap_id\n")

        # edge case 1: mix of True and False matches
        cl_match_ID_s = pd.Series(["depmap_id1", "depmap_id2", "error_depmap_id", "error_depmap_id", "depmap_id5"], index = ['ENSG1', 'ENSG2', 'ENSG3', "ENSG4", "ENSG5"])
        cl_match_ID_s.name = "top_corr_depmap_ID"
        cl_match_ID_s.index.name = "sample_id"
        logger.debug("\ncl_match_ID_s\n{}".format(cl_match_ID_s))
        logger.debug("\ncl_match_ID_s.index.name\n{}".format(cl_match_ID_s.index.name))

        query_expr_col_meta_df = pd.DataFrame({"sample_id":["ENSG1", "ENSG2", "ENSG3", "ENSG4", "ENSG5"], 
                                    "a":["1234", "5678", "2468", "3579", "2222"], 
                                    "b":["1234", "5678", "2468", "3579", "2222"], 
                                    "DepMap_ID":["depmap_id1", "depmap_id2", "depmap_id3", "depmap_id4", "depmap_id5"]})
        query_expr_col_meta_df = query_expr_col_meta_df.set_index("sample_id")
        logger.debug("\nquery_expr_col_meta_df.index.name:\n{}".format(query_expr_col_meta_df.index.name))

        logger.debug("\nquery_expr_col_meta_df\n{}".format(query_expr_col_meta_df))

        result_df_match = cla.compare_depmap_id(query_expr_col_meta_df, cl_match_ID_s)
        logger.debug("\nresult_df_match\n{}".format(result_df_match))

        expected_df_match = pd.DataFrame({"sample_id":["ENSG1", "ENSG2", "ENSG3", "ENSG4", "ENSG5"], 
                                    "a":["1234", "5678", "2468", "3579", "2222"], 
                                    "b":["1234", "5678", "2468", "3579", "2222"], 
                                    "DepMap_ID":["depmap_id1", "depmap_id2", "depmap_id3", "depmap_id4", "depmap_id5"],
                                    "top_corr_depmap_ID":["depmap_id1", "depmap_id2", "error_depmap_id", "error_depmap_id", "depmap_id5"],
                                    "cell_line_match":[True, True, False, False, True]})
        expected_df_match = expected_df_match.set_index("sample_id")
        logger.debug("\nexpected_df_match\n{}".format(expected_df_match))

        isEqual = expected_df_match == result_df_match
        logger.debug("\nisEqual:\n{}".format(isEqual))
        self.assertTrue(isEqual.all)

        # edge case 2: all False matches
        cl_match_ID_s = pd.Series(["error_depmap_id", "error_depmap_id", "error_depmap_id", "error_depmap_id", "error_depmap_id"], index = ['ENSG1', 'ENSG2', 'ENSG3', "ENSG4", "ENSG5"])
        cl_match_ID_s.name = "top_corr_depmap_ID"
        cl_match_ID_s.index.name = "sample_id"
        logger.debug("\ncl_match_ID_s\n{}".format(cl_match_ID_s))
        logger.debug("\ncl_match_ID_s.index.name\n{}".format(cl_match_ID_s.index.name))

        query_expr_col_meta_df = pd.DataFrame({"sample_id":["ENSG1", "ENSG2", "ENSG3", "ENSG4", "ENSG5"], 
                                    "a":["1234", "5678", "2468", "3579", "2222"], 
                                    "b":["1234", "5678", "2468", "3579", "2222"], 
                                    "DepMap_ID":["depmap_id1", "depmap_id2", "depmap_id3", "depmap_id4", "depmap_id5"]})
        query_expr_col_meta_df = query_expr_col_meta_df.set_index("sample_id")
        logger.debug("\nquery_expr_col_meta_df.index.name:\n{}".format(query_expr_col_meta_df.index.name))

        logger.debug("\nquery_expr_col_meta_df\n{}".format(query_expr_col_meta_df))

        result_df_match = cla.compare_depmap_id(query_expr_col_meta_df, cl_match_ID_s)
        logger.debug("\nresult_df_match\n{}".format(result_df_match))

        expected_df_match = pd.DataFrame({"sample_id":["ENSG1", "ENSG2", "ENSG3", "ENSG4", "ENSG5"], 
                                    "a":["1234", "5678", "2468", "3579", "2222"], 
                                    "b":["1234", "5678", "2468", "3579", "2222"], 
                                    "DepMap_ID":["depmap_id1", "depmap_id2", "depmap_id3", "depmap_id4", "depmap_id5"],
                                    "top_corr_depmap_ID":["error_depmap_id", "error_depmap_id", "error_depmap_id", "error_depmap_id", "error_depmap_id"],
                                    "cell_line_match":[False, False, False, False, False]})
        expected_df_match = expected_df_match.set_index("sample_id")
        logger.debug("\nexpected_df_match\n{}".format(expected_df_match))

        isEqual = expected_df_match == result_df_match
        logger.debug("\nisEqual:\n{}".format(isEqual))
        self.assertTrue(isEqual.all)        

        
    def test_check_for_warnings(self):
        logger.debug("\ntest_check_for_warnings\n")
        
        # edge case 1: no warnings/exceptions printed
        compared_depmaps_df = pd.DataFrame({"sample_id":["ENSG1", "ENSG2", "ENSG3", "ENSG4", "ENSG5"], 
                                    "a":["1234", "5678", "2468", "3579", "2222"], 
                                    "b":["1234", "5678", "2468", "3579", "2222"], 
                                    "DepMap_ID":["depmap_id1", "depmap_id2", "depmap_id3", "depmap_id4", "depmap_id5"],
                                    "top_corr_depmap_ID":["depmap_id1", "depmap_id2", "depmap_id3", "depmap_id4", "depmap_id5"],
                                    "cell_line_match":[True, True, True, True, True]})
        compared_depmaps_df = compared_depmaps_df.set_index("sample_id")
        cla.check_for_warnings(compared_depmaps_df)

        # edge case 2: mismatch; few warnings
        with self.assertRaises(cla.FhtbioinfpyCellLineAuthenticationNoDepMapIDMatch) as context:
            compared_depmaps_df = pd.DataFrame({"sample_id":["ENSG1", "ENSG2", "ENSG3", "ENSG4", "ENSG5"], 
                                    "a":["1234", "5678", "2468", "3579", "2222"], 
                                    "b":["1234", "5678", "2468", "3579", "2222"], 
                                    "DepMap_ID":["depmap_id1", "depmap_id2", "depmap_id3", "depmap_id4", "depmap_id5"],
                                    "top_corr_depmap_ID":["depmap_id1", "depmap_id2", "error_depmap_id", "error_depmap_id", "depmap_id5"],
                                    "cell_line_match":[True, True, False, False, True]})
            compared_depmaps_df = compared_depmaps_df.set_index("sample_id")
            cla.check_for_warnings(compared_depmaps_df)
                
        # edge case 3: all warnings
        with self.assertRaises(cla.FhtbioinfpyCellLineAuthenticationNoDepMapIDMatch) as context:
            compared_depmaps_df = pd.DataFrame({"sample_id":["ENSG1", "ENSG2", "ENSG3", "ENSG4", "ENSG5"], 
                                    "a":["1234", "5678", "2468", "3579", "2222"], 
                                    "b":["1234", "5678", "2468", "3579", "2222"], 
                                    "DepMap_ID":["depmap_id1", "depmap_id2", "depmap_id3", "depmap_id4", "depmap_id5"],
                                    "top_corr_depmap_ID":["error_depmap_id", "error_depmap_id", "error_depmap_id", "error_depmap_id", "error_depmap_id"],
                                    "cell_line_match":[False, False, False, False, False]})
            compared_depmaps_df = compared_depmaps_df.set_index("sample_id")
            cla.check_for_warnings(compared_depmaps_df)
           
           
    def test_save_compared_depmap_df(self):
            logger.debug("\ntest_save_df\n")

            #opening a file or directory --> make sure it gets cleaned up
            with tempfile.TemporaryDirectory(prefix = "fhtbioinfpy_test_cell_line_authentication") as tmpdirname:
                expected_output_file = cla.save_compared_depmap_df(self.df, self.exp_id, tmpdirname)
                logger.debug("\nexpected_output_file\n{}".format(expected_output_file))
                self.assertTrue(os.path.exists(expected_output_file))
                loaded_csv_file = pd.read_csv(expected_output_file, sep = "\t")
                self.assertFalse(loaded_csv_file.empty)
                self.assertEqual(loaded_csv_file.shape, (2, 4))
                self.assertEqual(loaded_csv_file.iloc[1,2], -1.0)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()



