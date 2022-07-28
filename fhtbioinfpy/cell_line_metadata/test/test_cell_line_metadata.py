from cgi import test
from re import A
import unittest
import logging
import fhtbioinfpy.setup_logger as setup_logger
import fhtbioinfpy.cell_line_metadata.cell_line_metadata as clm
import pandas as pd
import os
import tempfile
import numpy as np

logger = logging.getLogger(setup_logger.LOGGER_NAME)

class TestCellLineMetadata(unittest.TestCase):

    def setUp(self):
        logger.debug("setUp")   
        # In-house RNA-seq dataset
        self.query_expr_file = "./input_data/2022-02-10_heatmap_log2_TPM_r113153x1115.gctx"
        # Read the sample_info containing the DepMap_ID and the cell line name
        self.sample_info_file = "./input_data/sample_info.csv"
        self.shortened_sample_info_file = "./input_data/shortened_for_test_sample_info.csv"

    def tearDown(self):
        logger.debug("tearDown")


    def test_main_functional(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            logger.debug("test_main tmpdirname:  {}".format(tmpdirname))
            
            arg_list = ["--query_expr_file", self.query_expr_file,
                    "--sample_info_file", self.sample_info_file
            ]

            args = clm.build_parser().parse_args(arg_list)
            logger.debug(args)
            clm.main(args)
            expected_output_file = os.path.join(tmpdirname, "cell_line_metadata6x33.txt") # from build_output_file_name
            logger.debug("expected_output_file: {}".format(expected_output_file))
            self.assertTrue(os.path.isfile(expected_output_file))

            loaded_csv_file = pd.read_csv(expected_output_file, sep = "\t")
            self.assertFalse(loaded_csv_file.empty)    
        

    def test_remove_mouse_genes(self):
        test_dfq = pd.DataFrame({"gene":["ENSG1", "ENSG2", "ENSMUSG3", "ENSG4"], 
                                    "a":["1234", "5678", "2468", "3579"], 
                                    "b":["1234", "5678", "2468", "3579"], 
                                    "c":["1234", "5678", "2468", "3579"]})
        test_dfq = test_dfq.set_index("gene")
        logger.debug("test_dfq: \n{}".format(test_dfq))
        new_dfq = clm.remove_mouse_genes(test_dfq)
        logger.debug("new_dfq: \n{}".format(new_dfq))
        expected_dfq = pd.DataFrame({"gene":["ENSG1", "ENSG2", "ENSG4"], 
                                    "a":["1234", "5678", "3579"], 
                                    "b":["1234", "5678", "3579"], 
                                    "c":["1234", "5678", "3579"]})
        expected_dfq = expected_dfq.set_index("gene")
        logger.debug("expected_dfq: \n{}".format(expected_dfq))
        isEqual = new_dfq == expected_dfq
        logger.debug('\nisEqual:\n{}'.format(isEqual))
        self.assertTrue(isEqual.all)

    def test_drop_nas_query_df(self):
        test_dfq = pd.DataFrame({"gene":["ENSG1", "ENSG2", "ENSMUSG3", "ENSG4"], 
                                    "a":["1234", "5678", "2468", "3579"], 
                                    "b":["1234", "5678", "2468", "3579"], 
                                    "c":[np.nan, np.nan, np.nan, np.nan]})
        test_dfq = test_dfq.set_index("gene")
        logger.debug("test_dfq: \n{}".format(test_dfq))
        new_dfq, na_samples= clm.drop_nas_query_df(test_dfq)
        logger.debug("new_dfq: \n{}".format(new_dfq))
        expected_dfq = pd.DataFrame({"gene":["ENSG1", "ENSG2", "ENSMUSG3", "ENSG4"], 
                                    "a":["1234", "5678", "2468", "3579"], 
                                    "b":["1234", "5678", "2468", "3579"]})
        expected_dfq = expected_dfq.set_index("gene")
        logger.debug("expected_dfq: \n{}".format(expected_dfq))
        isEqual = new_dfq == expected_dfq
        logger.debug('\nisEqual:\n{}'.format(isEqual))
        self.assertTrue(isEqual.all)


    def test_read_metadata_col_gctx(self):
        test_gct_col = pd.DataFrame({"cid":["ENSG1", "ENSG2", "ENSMUSG3", "ENSG4"], 
                                    "bio_context_id":["1234", "5678", "2468", "3579"]})
        pass

    def test_drop_nas_metadata(self):
        test_dfq = pd.DataFrame({"gene":["ENSG1", "ENSG2", "ENSMUSG3", "ENSG4"], 
                                    "a":["1234", "5678", "2468", "3579"], 
                                    "b":["1234", "5678", "2468", "3579"], 
                                    "c":[np.nan, np.nan, np.nan, np.nan]})
        test_dfq = test_dfq.set_index("gene")
        logger.debug("test_dfq: \n{}".format(test_dfq))
        new_dfq= clm.drop_nas_metadata(test_dfq)
        logger.debug("new_dfq: \n{}".format(new_dfq))
        expected_dfq = pd.DataFrame({"gene":["ENSG1", "ENSG2", "ENSMUSG3", "ENSG4"], 
                                    "a":["1234", "5678", "2468", "3579"], 
                                    "b":["1234", "5678", "2468", "3579"],
                                    "c":["NOANNOTATION", "NOANNOTATION", "NOANNOTATION", "NOANNOTATION"]})
        expected_dfq = expected_dfq.set_index("gene")
        logger.debug("expected_dfq: \n{}".format(expected_dfq))
        isEqual = new_dfq == expected_dfq
        logger.debug('\nisEqual:\n{}'.format(isEqual))
        self.assertTrue(isEqual.all)

    '''    
    def test_add_cleaned_cl_name(self):
        test_dfq = pd.DataFrame({"gene":["ENSG1", "ENSG2", "ENSG3", "ENSG4"], 
                                    "a":["apple", "bana!na", "ORA,nge", "NCI.strawberry"], 
                                    "b":["1234", "5678", "2468", "3579"], 
                                    "c":["1234", "5678", "2468", "3579"]})
        test_dfq = test_dfq.set_index("gene")
        col_name = 'a'
        cleaned_dfq = clm.add_cleaned_cl_name(test_dfq, col_name)
        logger.debug('\ncleaned_df:\n{}'.format(cleaned_dfq))
        expected_dfq = pd.DataFrame({"gene":["ENSG1", "ENSG2", "ENSG3", "ENSG4"],
                                "a":["apple", "bana!na", "ORA,nge", "NCI.strawberry"],  
                                "b":["1234", "5678", "2468", "3579"], 
                                "c":["1234", "5678", "2468", "3579"],
                                "cleaned_a":["APPLE", "BANANA", "ORANGE", "STRAWBERRY"]})
        expected_dfq = expected_dfq.set_index("gene")
        logger.debug('\nexpected_df:\n{}'.format(expected_dfq))
        isEqual = cleaned_dfq == expected_dfq
        logger.debug('\nisEqual:\n{}'.format(isEqual))
        self.assertTrue(isEqual.all)
    '''
    '''
    def test_shorten_sample_info(self):
        info_path = self.shortened_sample_info_file
        new_sample_info = clm.shorten_sample_info(info_path)
        logger.debug("\nnew_sample_info:\n{}".format(new_sample_info))
        expected_sample_info = pd.DataFrame({"DepMap_ID":["ACH-000001", "ACH-000002", "ACH-000003", "ACH-000004", "ACH-000005"],
                              "stripped_cell_line_name": ["NIHOVCAR3", "HL60", "CACO2", "HEL", "HEL9217"]})
        logger.debug('\nexpected_sample_info:\n{}'.format(expected_sample_info))
        isEqual = new_sample_info == expected_sample_info
        logger.debug('\nisEqual:\n{}'.format(isEqual))
        self.assertTrue(isEqual.all)
    '''    

    def test_drop_nas_sample_info(self):
        #edge case 1: no NAs to be dropped in the sample info
        test_sample_info = pd.DataFrame({"DepMap_ID":["1234", "5678", "2468", "3579"], 
                           "stripped_cell_line_name":["ab1234", "ab5678", "ab2468", "ab3579"]})
        new_sample_info = clm.drop_nas_sample_info(test_sample_info)
        logger.debug('\nnew_sample_info:\n{}'.format(new_sample_info))
        expected_sample_info = pd.DataFrame({"DepMap_ID":["1234", "5678", "2468", "3579"], 
                               "stripped_cell_line_name":["ab1234", "ab5678", "ab2468", "ab3579"]})
        logger.debug('\nexpected_sample_info:\n{}'.format(expected_sample_info))
        isEqual = new_sample_info == expected_sample_info
        logger.debug('\nisEqual:\n{}'.format(isEqual))
        self.assertTrue(isEqual.all)

        #edge case 2: NAs in sample info are dropped
        test_sample_info = pd.DataFrame({"DepMap_ID":["1234", "5678", "2468", "3579"], 
                           "stripped_cell_line_name":[np.nan, "ab5678", "ab2468", np.nan]})
        logger.debug('\ntest_sample_info:\n{}'.format(test_sample_info))
        new_sample_info = clm.drop_nas_sample_info(test_sample_info)
        new_sample_info = new_sample_info.reset_index(drop=True)
        logger.debug('\nnew_sample_info:\n{}'.format(new_sample_info))        
        expected_sample_info = pd.DataFrame({"DepMap_ID":["5678", "2468"], 
                               "stripped_cell_line_name":["ab5678", "ab2468"]})
        logger.debug('\nexpected_sample_info:\n{}'.format(expected_sample_info))
        isEqual = new_sample_info == expected_sample_info
        logger.debug('\nisEqual:\n{}'.format(isEqual))
        self.assertTrue(isEqual.all)

    def test_merge_metadata_w_sample_info(self):
        # edge case 1: merges properly, metadata and sample_info columns are equal
        metadata = pd.DataFrame({"sample":["1234", "5678", "2468", "3579"], 
                         "bio_context_id":["1234", "5678", "2468", "3579"],
                 "cleaned_bio_context_id":["AB1234", "AB5678", "AB2468", "AB3579"]})
        logger.debug('\n metadata:\n{}'.format(metadata))

        sample_info = pd.DataFrame({"DepMap_ID":["1234", "5678", "2468", "3579"], 
                      "stripped_cell_line_name":["ab1234", "ab5678", "ab2468", "ab3579"],
              "cleaned_stripped_cell_line_name":["AB1234", "AB5678", "AB2468", "AB3579"]})
        logger.debug('\n sample_info:\n{}'.format(sample_info))

        merging_result = clm.merge_metadata_w_sample_info(metadata, sample_info)
        logger.debug('\nmerging_result:\n{}'.format(merging_result))
        expected_merging = pd.DataFrame({"sample":["1234", "5678", "2468", "3579"], 
                                 "bio_context_id":["1234", "5678", "2468", "3579"],
                         "cleaned_bio_context_id":["AB1234", "AB5678", "AB2468", "AB3579"],
                                      "DepMap_ID":["1234", "5678", "2468", "3579"], 
                        "stripped_cell_line_name":["ab1234", "ab5678", "ab2468", "ab3579"]})
        logger.debug('\nexpected_merging:\n{}'.format(expected_merging))
        isEqual = merging_result == expected_merging
        logger.debug('\nisEqual:\n{}'.format(isEqual))
        self.assertTrue(isEqual.all)

        #edge case 2: metadata column and sample_info column are not equal
        with self.assertRaises(clm.FhtbioinfpyCellLineMetadataCheckEquivalenceOfBCIDandStrippedCLName) as context:
            metadata = pd.DataFrame({"sample":["1234", "5678", "2468", "3579"], 
                            "bio_context_id":["1234", "5678", "2468", "3579"],
                    "cleaned_bio_context_id":["AB1234", "AB5678", "AB2468", "AB3579"]})
            logger.debug('\n metadata:\n{}'.format(metadata))

            sample_info = pd.DataFrame({"DepMap_ID":["1234", "5678", "2468", "3579"], 
                        "stripped_cell_line_name":["ab1234", "ab5678", "ab2468", "ab3579"],
                "cleaned_stripped_cell_line_name":["AB1234", "AB5678", "AB2468", "CD3579"]})
            logger.debug('\n sample_info:\n{}'.format(sample_info))
            clm.merge_metadata_w_sample_info(metadata, sample_info)

        
    def test_build_output_file_name(self):
        test_df = pd.DataFrame({"a":["1A", "2B", "CC", "3D"], "b":["5F", "GG", "4H", "HI"]})

        output_filename = clm.build_output_file_name(test_df.shape)
        logger.debug('\nOUTPUT FILENAME:\n{}'.format(output_filename))
        expected_filename = "cell_line_metadata4x2.txt"
        self.assertEqual(output_filename, expected_filename)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()