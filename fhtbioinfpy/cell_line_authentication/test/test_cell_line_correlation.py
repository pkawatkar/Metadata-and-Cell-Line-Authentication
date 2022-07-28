from cgi import test
from re import A
import unittest
import logging
import fhtbioinfpy.setup_logger as setup_logger
import fhtbioinfpy.cell_line_authentication.cell_line_correlation as clc
import pandas as pd
import os
import tempfile
import numpy as np
import cmapPy
import cmapPy.pandasGEXpress.parse as parse

logger = logging.getLogger(setup_logger.LOGGER_NAME)

class TestCellLineCorrelation(unittest.TestCase):

    def setUp(self):
        logger.debug("setUp")   
        
        self.query_expr_filepath = "./assets/test_cell_line_auth/test_cell_line_auth_query_data_r110x86.gctx"
        self.ref_expr_filepath   = "./assets/test_cell_line_auth/test_cell_line_auth_ref_data_r108x1406.gctx"

        self.df = pd.DataFrame({"c":[-1.0, 1.0], "d":[1.0, -1.0], "e":[1.0, -1.0]})
        self.df = self.df.set_index(pd.Series(['a', 'b']))

        self.query_expr_gctoo = parse.parse(self.query_expr_filepath)
        logger.debug('\nself.df.shape\n{}'.format(self.df.shape))
        logger.debug('\nself.df.head\n{}'.format(self.df.head()))
        
        self.exp_id = "test_experiment_id"
        self.subdir = "./assets/cell_line_auth/"
        self.row_metadata_for_matching = "gene_symbol"
        

    def tearDown(self):
        logger.debug("tearDown")


    def test_main_functional(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            logger.debug("test_main tmpdirname:  {}".format(tmpdirname))
            
            arg_list = ["--query_expr_filepath", self.query_expr_filepath,
                        "--ref_expr_filepath", self.ref_expr_filepath,
                        "--experiment_id", self.exp_id,
                        "--output_subdir", tmpdirname,
                        "--row_metadata_for_matching", self.row_metadata_for_matching
            ]

            args = clc.build_parser().parse_args(arg_list)
            logger.debug(args)
            clc.main(args)

            expected_output_file = os.path.join(tmpdirname, "test_experiment_id_cell_line_authentication_corr_r86x1406.txt")
            logger.debug("\nexpected_output_file: \n{}".format(expected_output_file))
            self.assertTrue(os.path.isfile(expected_output_file))
            loaded_csv_file = pd.read_csv(expected_output_file, sep = "\t")
            self.assertFalse(loaded_csv_file.empty)  

    def test_read_RNA_seq_gctx(self):
        logger.debug("\ntest_read_RNA_seq_gctx\n")

        self.assertTrue(os.path.exists(self.query_expr_filepath))
        # check the two dataframes are equal
        query_expr_gctoo = clc.read_RNA_seq_gctx(self.query_expr_filepath)
        query_expr_df = query_expr_gctoo.data_df

        logger.debug("\nquery_expr_data.shape: {}".format(query_expr_df.shape))   
        self.assertFalse(query_expr_df.empty)
        self.assertEqual(query_expr_df.shape, (110, 86)) 
        logger.debug("\nquery_expr_df.iloc[1,2]:\n{}".format(query_expr_df.iloc[1,2])) 
        self.assertEqual(query_expr_df.iloc[1,2], 3.307300090789795) 


    def test_load_ref_expr_data(self):
        logger.debug("\ntest_load_ref_expr_data\n")

        self.assertTrue(os.path.exists(self.ref_expr_filepath))
        # check the two dataframes are equal
        ref_expr_gctoo = clc.load_ref_expr_data(self.ref_expr_filepath)
        ref_expr_df = ref_expr_gctoo.data_df
        logger.debug("\nref_expr_df.shape: \n{}".format(ref_expr_df.shape))
        self.assertFalse(ref_expr_df.empty)
        self.assertEqual(ref_expr_df.shape, (108, 1406)) 
        logger.debug("\nref_expr_df.iloc[1,2]:\n{}".format(ref_expr_df.iloc[1,2])) 
        self.assertEqual(ref_expr_df.iloc[1,2], 4.372951984405518) 

    def test_query_set_index(self):
        logger.debug("\ntest_query_set_index\n")

        query_expr_df = clc.query_set_index(self.query_expr_gctoo, self.row_metadata_for_matching)
        #check gene_index spots, first spot, last spot
        self.assertTrue(query_expr_df.index.name==self.row_metadata_for_matching)
        self.assertFalse(query_expr_df.empty)
        logger.debug("\nquery_expr_df.shape: \n{}".format(query_expr_df.shape))
        self.assertEqual(query_expr_df.shape, (110, 86)) 
        logger.debug("\nquery_expr_df.iloc[1,2]:\n{}".format(query_expr_df.iloc[1,2])) 
        self.assertEqual(query_expr_df.iloc[1,2], 3.307300090789795) 

    def test_build_matched_datasets(self):
        logger.debug("\ntest_build_matched_datasets\n")

        # edge case 1 
        test_query_expr = pd.DataFrame({"gene":["ENSG1", "ENSG2", "ENSG3", "ENSG4"], 
                                    "a":["1234", "5678", "2468", "3579"], 
                                    "b":["1234", "5678", "2468", "3579"], 
                                    "c":["1234", "5678", "2468", "3579"]})
        test_query_expr = test_query_expr.set_index("gene")
        test_ref_expr = pd.DataFrame({"gene":["ENSG1", "ENSG2", "ENSMUSG3", "ENSG4"], 
                                    "a":["1234", "5678", "2468", "3579"], 
                                    "b":["1234", "5678", "2468", "3579"], 
                                    "c":["1234", "5678", "2468", "3579"]})
        test_ref_expr = test_ref_expr.set_index("gene")

        expected_qe = pd.DataFrame({"gene":["ENSG1", "ENSG2", "ENSG4"], 
                                    "a":["1234", "5678", "3579"], 
                                    "b":["1234", "5678", "3579"], 
                                    "c":["1234", "5678", "3579"]})
        expected_qe = expected_qe.set_index("gene")
        expected_re = pd.DataFrame({"gene":["ENSG1", "ENSG2", "ENSG4"], 
                                    "a":["1234", "5678", "3579"], 
                                    "b":["1234", "5678", "3579"], 
                                    "c":["1234", "5678", "3579"]})
        expected_re = expected_re.set_index("gene")

        new_query_expr, new_ref_expr = clc.build_matched_datasets(test_query_expr, test_ref_expr)

        query_expEqual = expected_qe == new_query_expr
        ref_expEqual   = expected_re == new_ref_expr

        self.assertTrue(query_expEqual.all)
        self.assertTrue(ref_expEqual.all)

        # edge case 2 - different dataframe sizes
        test_query_expr = pd.DataFrame({"gene":["ENSG1", "ENSG2", "ENSG3", "ENSG4", "ENSG5"], 
                                    "a":["1234", "5678", "2468", "3579", "2222"], 
                                    "b":["1234", "5678", "2468", "3579", "2222"], 
                                    "c":["1234", "5678", "2468", "3579", "2222"]})
        test_query_expr = test_query_expr.set_index("gene")
        test_ref_expr = pd.DataFrame({"gene":["ENSG1", "ENSG2", "ENSMUSG3", "ENSG4"], 
                                    "a":["1234", "5678", "2468", "3579"], 
                                    "b":["1234", "5678", "2468", "3579"], 
                                    "c":["1234", "5678", "2468", "3579"]})
        test_ref_expr = test_ref_expr.set_index("gene")

        expected_qe = pd.DataFrame({"gene":["ENSG1", "ENSG2", "ENSG4"], 
                                    "a":["1234", "5678", "3579"], 
                                    "b":["1234", "5678", "3579"], 
                                    "c":["1234", "5678", "3579"]})
        expected_qe = expected_qe.set_index("gene")
        expected_re = pd.DataFrame({"gene":["ENSG1", "ENSG2", "ENSG4"], 
                                    "a":["1234", "5678", "3579"], 
                                    "b":["1234", "5678", "3579"], 
                                    "c":["1234", "5678", "3579"]})
        expected_re = expected_re.set_index("gene")

        new_query_expr, new_ref_expr = clc.build_matched_datasets(test_query_expr, test_ref_expr)

        query_expEqual = expected_qe == new_query_expr
        ref_expEqual = expected_re == new_ref_expr

        self.assertTrue(query_expEqual.all)
        self.assertTrue(ref_expEqual.all)
        

    def test_run_correlation_calculation(self):
        logger.debug("\ntest_run_correlation_calculation\n")

        matched_query_df = pd.DataFrame({"a": [2, 3], "b": [7, 5]})
        matched_ref_df = pd.DataFrame({"c": [13, 11], "d": [17, 19], "e": [23, 29]})
        result = clc.run_correlation_calculation(matched_query_df, matched_ref_df)
        logger.debug("\nfast corr\n {}".format(result))

        expected = pd.DataFrame({"c":[-1.0, 1.0], "d":[1.0, -1.0], "e":[1.0, -1.0]})
        expected = expected.set_index(pd.Series(['a', 'b']))
        logger.debug("\nexpected\n{}".format(expected))

        #check size/shape:
        result_shape = result.shape
        logger.debug("\nresult_shape: \n{}".format(result_shape))
        expected_shape = expected.shape
        logger.debug("\nexpected_shape: \n{}".format(expected_shape))
        self.assertEqual(result_shape, expected_shape)

        #check values of dataframe:
        valuesEqual = result == expected
        self.assertTrue(valuesEqual.all)
        

    def test_save_corr_df(self):
        logger.debug("\ntest_save_df\n")

        #opening a file or directory --> make sure it gets cleaned up
        with tempfile.TemporaryDirectory(prefix = "fhtbioinfpy_test_cell_line_authentication") as tmpdirname:
            expected_output_file = clc.save_corr_df(self.df, self.exp_id, tmpdirname)
            logger.debug("\nexpected_output_file\n{}".format(expected_output_file))
            self.assertTrue(os.path.exists(expected_output_file))

            loaded_csv_file = pd.read_csv(expected_output_file, sep = "\t")
            self.assertFalse(loaded_csv_file.empty)
            self.assertEqual(loaded_csv_file.shape, (2, 4))
            self.assertEqual(loaded_csv_file.iloc[1,2], -1.0)
            

if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
