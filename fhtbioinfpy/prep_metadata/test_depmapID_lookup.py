
import unittest
import logging
import fhtbioinfpy.setup_logger as setup_logger
import fhtbioinfpy.prep_metadata.depmapID_lookup as dmid
import pandas as pd
import os
import tempfile



logger = logging.getLogger(setup_logger.LOGGER_NAME)

class TestPrepMetadata(unittest.TestCase):

    def setUp(self):
        logger.debug("Setup")
        
        self.sample_info_file = "./input_data/sample_info.csv"
        self.shortened_sample_info_file = "./input_data/shortened_for_test_sample_info.csv"

    def tearDown(self):
        logger.debug("Teardown")
        
    def test_add_cleaned_cl_name(self):
        test_dfq = pd.DataFrame({"gene":["ENSG1", "ENSG2", "ENSG3", "ENSG4"], 
                                    "a":["apple", "bana!na", "ORA,nge", "NCI.strawberry"], 
                                    "b":["1234", "5678", "2468", "3579"], 
                                    "c":["1234", "5678", "2468", "3579"]})
        test_dfq = test_dfq.set_index("gene")
        col_name = 'a'
        cleaned_dfq = dmid.add_cleaned_cl_name(test_dfq, col_name)
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

    def test_shorten_sample_info(self):
        info_path = self.shortened_sample_info_file
        new_sample_info = dmid.shorten_sample_info(info_path)
        logger.debug("\nnew_sample_info:\n{}".format(new_sample_info))
        expected_sample_info = pd.DataFrame({"DepMap_ID":["ACH-000001", "ACH-000002", "ACH-000003", "ACH-000004", "ACH-000005"],
                              "stripped_cell_line_name": ["NIHOVCAR3", "HL60", "CACO2", "HEL", "HEL9217"]})
        logger.debug('\nexpected_sample_info:\n{}'.format(expected_sample_info))
        isEqual = new_sample_info == expected_sample_info
        logger.debug('\nisEqual:\n{}'.format(isEqual))
        self.assertTrue(isEqual.all)


    def test_verify_match(self):
        #edge case 1: all stripped_bio_context_id values and stripped_cell_line_name values match and corresponding depmapIDs are found
        metadata = pd.DataFrame({"a":["1234", "5678", "2468", "3579"], 
                                    "b":["1234", "5678", "2468", "3579"], 
                                    "c":["1234", "5678", "2468", "3579"],
                                    "stripped_bio_context_id":["ENSG1", "ENSG2", "ENSG3", "ENSG4"]})
        sample_info = pd.DataFrame({"stripped_cell_line_name":["ENSG1", "ENSG2", "AB78", "AB80", "AB98", "AB09", "ENSG3", "ENSG4"], 
                                    "DepMap_ID":["1234", "5678", "2468", "3579", "3456", "4567", "6789", "9876"]})
        
        result = dmid.verify_match(metadata, sample_info)
        expected_result = pd.DataFrame({"a":["1234", "5678", "2468", "3579"], 
                                    "b":["1234", "5678", "2468", "3579"], 
                                    "c":["1234", "5678", "2468", "3579"],
                                    "stripped_bio_context_id":["ENSG1", "ENSG2", "ENSG3", "ENSG4"],
                                    "DepMap_ID":["1234", "5678", "6789", "9876"]})
        isEqual = result == expected_result
        self.assertTrue(isEqual.all)


        #edge case 2:  one bio_context_id values that don't match and have a corresponding DepMap_ID
        with self.assertRaises(dmid.FhtbioinfpydepmapIDVerifyEquivalenceOfCLMandBCID) as context:
            metadata = pd.DataFrame({"a":["1234", "5678", "2468", "3579"], 
                                    "b":["1234", "5678", "2468", "3579"], 
                                    "c":["1234", "5678", "2468", "3579"],
                                    "stripped_bio_context_id":["ENSG1", "ENSG2", "ENSG3", "ENSG4"]})
            sample_info = pd.DataFrame({"stripped_cell_line_name":["ENSG1", "ENSG5", "AB78", "AB80", "AB98", "AB09", "ENSG3", "ENSG4"], 
                                        "DepMap_ID":["1234", "5678", "2468", "3579", "3456", "4567", "6789", "9876"]})
            
            dmid.verify_match(metadata, sample_info)

        

         #edge case 3: more than one bio_context_id values that don't match and have a corresponding DepMap_ID 
        with self.assertRaises(dmid.FhtbioinfpydepmapIDVerifyEquivalenceOfCLMandBCID) as context:
            metadata = pd.DataFrame({"a":["1234", "5678", "2468", "3579"], 
                                    "b":["1234", "5678", "2468", "3579"], 
                                    "c":["1234", "5678", "2468", "3579"],
                                    "stripped_bio_context_id":["ENSG10", "ENSG20", "ENSG30", "ENSG40"]})
            sample_info = pd.DataFrame({"stripped_cell_line_name":["ENSG1", "ENSG5", "AB78", "AB80", "AB98", "AB09", "ENSG3", "ENSG4"], 
                                        "DepMap_ID":["1234", "5678", "2468", "3579", "3456", "4567", "6789", "9876"]})
            
            dmid.verify_match(metadata, sample_info)

 



if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()