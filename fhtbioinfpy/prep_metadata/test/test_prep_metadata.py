import unittest
import logging
import fhtbioinfpy.setup_logger as setup_logger
import fhtbioinfpy.prep_metadata.prep_metadata as pm
import pandas as pd
import os
import tempfile
import fhtbioinfpy.prep_metadata.depmapID_lookup as dmid



logger = logging.getLogger(setup_logger.LOGGER_NAME)

class TestPrepMetadata(unittest.TestCase):

    def setUp(self):
        logger.debug("Setup")
        self.metadata_subdir = "./assets/metadata/"
        self.metadata_file ="2021-11-09-next-gen-sequencing-annotations-ATACseq-KA.txt"
        self.input_metadata_file = "./assets/metadata/2021-11-09-next-gen-sequencing-annotations-ATACseq-KA.txt"
        
        #self.metadata_file ="2022-04-13-NS-22.0019-metadata-OP-KYSE70-shRNA-TP63-knockdowns.txt"
        #self.input_metadata_file = "./assets/metadata/2022-04-13-NS-22.0019-metadata-OP-KYSE70-shRNA-TP63-knockdowns.txt"
        


        #2022-04-13-NS-22.0019-metadata-OP-KYSE70-shRNA-TP63-knockdowns
        self.df = pd.read_csv(self.input_metadata_file, sep="\t", index_col = "sample_id")
        logger.debug('self.df.shape{}'.format(self.df.shape))
        logger.debug('self.df.head{}'.format(self.df.head()))
        # encoding="ISO-8859-1"

    def tearDown(self):
        logger.debug("Teardown")

    def test_main_functional(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            logger.debug("test_main tmpdirname:  {}".format(tmpdirname))
            arg_list = ["--input_metadata_file", self.input_metadata_file,
                        "--experiment_id", "test_main_experiment_id",
                        "--metadata_columns_to_build_groups", "pert_id", "bio_context_id",
                        "--output_metadata_subdir", tmpdirname,
                        "--sample_info_file",  "./input_data/sample_info.csv"
            ]

            args = pm.build_parser().parse_args(arg_list)
            logger.debug(args)
            pm.main(args)

            expected_output_file = os.path.join(tmpdirname, "test_main_experiment_id_metadata_r6x33.txt") # from build_output_file_name
            logger.debug("expected_output_file: {}".format(expected_output_file))
            self.assertTrue(os.path.isfile(expected_output_file))

            loaded_csv_file = pd.read_csv(expected_output_file, sep = "\t")
            self.assertFalse(loaded_csv_file.empty)

    def test_load_metadata(self):
        logger.debug("test_load_metadata")
        # check that the input_metadata_file exists
        self.assertTrue(os.path.exists(self.input_metadata_file))
        # check the two dataframes are equal
        orig_metadata = pm.load_metadata(self.input_metadata_file)
        logger.debug("orig_metadata.shape: {}".format(orig_metadata.shape))
        # self.assertEqual(orig_metadata.shape, self.df.shape)

    def test_remove_samples(self):
        logger.debug("test_remove_samples")
        #edge case 1 --> remove 1 row
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "c":["happy","sad", "laugh"], "d":["python","java", "html"]})
        logger.debug('\ntest_df:\n{}'.format(test_df))
        rows_to_remove = 1
        updated_df = pm.remove_samples(test_df, rows_to_remove)
        logger.debug('\nupdated_df:\n{}'.format(updated_df))
        expected_df = pd.DataFrame({"a":["hello", "basketball"], "b":["yellow", "book"], "c":["happy", "laugh"], "d":["python", "html"]})
        logger.debug('\nexpected_df:\n{}'.format(expected_df))
 
        isEqual = updated_df == expected_df
        logger.debug('\nisEqual:\n{}'.format(isEqual))
        self.assertTrue(isEqual.all)
    
        #edge case 2 --> remove multiple rpws
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball", "banana"], "b":["yellow","cherry", "book", "carrot"], "c":["happy","sad", "laugh", "angry"], "d":["python","java", "html", "rstudio"]})
        logger.debug('\ntest_df:\n{}'.format(test_df))
        rows_to_remove = [1, 3]
        updated_df = pm.remove_samples(test_df, rows_to_remove)
        logger.debug('\nupdated_df:\n{}'.format(updated_df))
        expected_df = pd.DataFrame({"a":["hello", "basketball"], "b":["yellow", "book"], "c":["happy", "laugh"], "d":["python", "html"]})
        logger.debug('\nexpected_df:\n{}'.format(expected_df))
        
        isEqual = updated_df == expected_df
        logger.debug('\nisEqual:\n{}'.format(isEqual))
        self.assertTrue(isEqual.all)
        
    def test_check_replicate_number_per_group(self):
        logger.debug("test_check_replicate_number_per_group")
        #edge case 1 --> 3 replicates
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball", "apple", "crayon", "dog", "apple", "crayon", "dog"], "b":["yellow","cherry", "book", "yellow","cherry", "book", "yellow","cherry", "book"], "c":["happy","sad", "laugh", "apple", "crayon", "dog", "apple", "crayon", "dog"]})
        selected_cols = "b"
        all3s = pm.check_replicate_number_per_group(test_df, selected_cols, 3)
        logger.debug('\nall3s:\n{}'.format(all3s))
        checker = all3s == 3
        logger.debug('\nchecker:\n{}'.format(checker))
        isequal = checker.all()
        logger.debug('\nisequal:\n{}'.format(isequal))
        self.assertTrue(isequal.all())

        #edge case 2 --> replicates are not what were expected and user chooses to print exception error and break code
        with self.assertRaises(pm.FhtbioinfpyPrepMetadataCheckValueReplicatesForGroups) as context:
            test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "c":["happy","sad", "laugh"]})
            selected_cols = "b"
            pm.check_replicate_number_per_group(test_df, selected_cols, 3, lambda : True)

        #edge case 3 --> replicates are not what were expected but user ignores exception error and continues
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["sample1","sample2", "sample3"], "c":["happy","sad", "laugh"]})
        selected_cols = "b"
        try:
            pm.check_replicate_number_per_group(test_df, selected_cols, 3, lambda : False)
        except pm.FhtbioinfpyPrepMetadataCheckValueReplicatesForGroups:
            self.fail("check_replicate_number_per_group() unexpectedly raised FhtbioinfpyPrepMetadataCheckValueReplicatesForGroups exception!")

    def test_add_string_concated_columns(self):
        logger.debug("test_add_string_concated_columns")
        test_df = pd.DataFrame({"pert_dose":[3, 4, 5], "pert_dose_unit": ["nM", "nM", "nM"], "pert_time":[10, 11, 12], "pert_time_unit": ["h", "h", "h"]})
        logger.debug('\ntest_df:\n{}'.format(test_df))
        new_df = pm.add_string_concated_columns(test_df)
        logger.debug('\nnew_df:\n{}'.format(new_df))
        subset = False
        if set(['pert_dose_str','pert_time_str']).issubset(new_df.columns):
            subset = True
        logger.debug('\nsubset:\n{}'.format(subset))
        self.assertTrue(subset)
        self.assertEqual(new_df["pert_dose_str"][0],"3 nM")

    def test_verify_group_def_columns_in_metadata(self):
        logger.debug("test_verify_group_def_columns_in_metadata")

        #edge case 1 --> only one column in df is selected
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "c":["happy","sad", "laugh"], "d":["python","java", "html"]})
        col = "c"
        found = pm.verify_group_def_columns_in_metadata(test_df, col)
        logger.debug('\nfound:\n{}'.format(found))
        self.assertTrue(found)

        #edge case 2 --> more than one column in df is selected
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "c":["happy","sad", "laugh"], "d":["python","java", "html"]})
        cols = ["a", "d"]
        found = pm.verify_group_def_columns_in_metadata(test_df, cols)
        logger.debug('\nfound:\n{}'.format(found))
        self.assertTrue(found)

        with self.assertRaises(pm.FhtbioinfpyPrepMetadataVerifyColumnsForGroups) as context:
            #edge case 3 --> one column in df is selected but not found in df (print out exception)
            test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "c":["happy","sad", "laugh"], "d":["python","java", "html"]})
            cols = ["a", "f"] #f is a column not found in df
            pm.verify_group_def_columns_in_metadata(test_df, cols)
        with self.assertRaises(pm.FhtbioinfpyPrepMetadataVerifyColumnsForGroups) as context:
            #edge case 4 --> more than one column not found in df (print out exception)
            test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "c":["happy","sad", "laugh"], "d":["python","java", "html"]})
            cols = ["a", "f", "h"] #columns f and h are not found in df
            pm.verify_group_def_columns_in_metadata(test_df, cols)  


    def test_build_group_names(self):
        logger.debug("test_build_group_names")
        #edge case 1 --> one selected column
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "c":["happy","sad", "laugh"], "d":["python","java", "html"]})
        logger.debug('\ntest_df:\n{}'.format(test_df))
        selected_cols = ["c"]
        new_groups_df = pm.build_group_names(test_df, selected_cols)
        logger.debug('\nnew_groups_df:\n{}'.format(new_groups_df))
        expected_group_names = ["happy", "sad", "laugh"]
        logger.debug('\nexpected_group_names:\n{}'.format(expected_group_names))
        self.assertEqual(list(new_groups_df), expected_group_names)

        #edge case 2 --> more than one selected column
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "c":["happy","sad", "laugh"], "d":["python","java", "html"]})
        logger.debug('\ntest_df:\n{}'.format(test_df))
        selected_cols = ["a","c"]
        new_groups_df = pm.build_group_names(test_df, selected_cols)
        logger.debug('\nnew_groups_df:\n{}'.format(list(new_groups_df)))
        expected_group_names = ["hello_happy", "apple_sad", "basketball_laugh"]
        logger.debug('\nexpected_groups_df:\n{}'.format(expected_group_names))
        self.assertEqual(list(new_groups_df), expected_group_names)
    
    def test_add_groups_to_df(self):
        logger.debug("test_add_groups_to_df")
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"]})
        logger.debug('\ntest_df:\n{}'.format(test_df))
        groups_df = pd.DataFrame({"group":["sample1_data", "sample2_data", "sample3_data"]})
        logger.debug('\ngroups_df:\n{}'.format(groups_df))
        expected_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "group":["sample1_data", "sample2_data", "sample3_data"]})
        logger.debug('\nexpected_df:\n{}'.format(expected_df))
        updated_df = pm.add_groups_to_df(test_df, groups_df)
        logger.debug('\nupdated_df:\n{}'.format(updated_df))
        isEqual = updated_df == expected_df
        logger.debug('\nisEqual:\n{}'.format(isEqual))
        final_bool = isEqual.all
        logger.debug('\nfinal_bool:\n{}'.format(final_bool))
        self.assertTrue(final_bool)

    def test_add_experiment_id(self):
        logger.debug("test_add_experiment_id")
        experiment_id = "NS-21.0138"
        new_df = pm.add_experiment_id(self.df, experiment_id)
        subset = {'experiment_id'}.issubset(new_df)
        self.assertTrue(subset)

    def test_create_R_Groups(self):
        logger.debug("test_create_R_Groups")
        
        #edge case 1 --> all start with a number / all are True
        test_df = pd.DataFrame({"a":[str(3), str(4), str(5)], "b":[str(8), str(9), str(10)]})
        logger.debug('\ntest_df:\n{}'.format(test_df))
        new_df = pm.create_R_Groups(test_df, "a")
        logger.debug('\nnew_df:\n{}'.format(new_df))
        check1 = new_df['R_group'].str[0:2]#first two characters of the string in column 'a'
        logger.debug('\ncheck1:\n{}'.format(check1))
        check2 = check1 == "x_"
        logger.debug('\ncheck2:\n{}'.format(check2))
        self.assertTrue(check2.all()) #all() true if all true

        # edge case 2 --> none start with a number / all are False
        test_df = pd.DataFrame({"a":["AB", "BC", "CD", "DE"], "b":["EF", "FG", "GH", "HI"]})
        logger.debug('test_df: \n{}\n'.format(test_df))
        new_df = pm.create_R_Groups(test_df, "b")
        logger.debug('new_Df: \n{}\n'.format(new_df))
        check1 = new_df.b.str[0:2]
        logger.debug('\ncheck1:\n{}'.format(check1))
        check2 = check1 == "x_"
        logger.debug('\ncheck2:\n{}'.format(check2))
        self.assertFalse(check2.all())
        
        # edge case 3 --> mix of strings that start with numbers and letters
        test_df = pd.DataFrame({"a":["1A", "2B", "CC", "3D"], "b":["5F", "GG", "4H", "HI"]})
        logger.debug('\ntest_df:\n{}'.format(test_df))
        new_df = pm.create_R_Groups(test_df, "b")
        allWork = True
        check1 = new_df['R_group'].str[0:2]
        logger.debug('\ncheck1:\n{}'.format(check1))
        check2 = check1 =="x_"
        logger.debug('\ncheck2:\n{}'.format(check2))
        if check1[0]=='^[0-9]': #after the method is called none of the strings should start with a number
            allWork = False
            self.assertEquals(check1, "x_") #should not pass if this if statement is run
            
        #else--> check1 is a nonnumber string and isn't equal to "x_" so continues running through rows
        self.assertTrue(allWork)
    
    def test_build_output_file_name(self):
        experiment_id = "id"
        test_df = pd.DataFrame({"a":["1A", "2B", "CC", "3D"], "b":["5F", "GG", "4H", "HI"]})

        output_filename = pm.build_output_file_name(experiment_id, test_df.shape)
        logger.debug('\nOUTPUT FILENAME:\n{}'.format(output_filename))
        expected_filename = "id_metadata_r4x2.txt"
        self.assertEqual(output_filename, expected_filename)

        


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
