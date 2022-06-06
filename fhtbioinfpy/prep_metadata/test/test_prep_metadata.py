import unittest
import logging
import fhtbioinfpy.setup_logger as setup_logger
import fhtbioinfpy.prep_metadata.prep_metadata as pm
import pandas as pd
import os.path
import tempfile


logger = logging.getLogger(setup_logger.LOGGER_NAME)

class TestPrepMetadata(unittest.TestCase):

    def setUp(self):
        print("Setup")
        self.metadata_subdir = "./assets/metadata/"
        self.metadata_file ="2021-11-09 next-gen-sequencing-annotations-ATACseq-KA.txt"
        self.input_metadata_file = "./assets/metadata/2021-11-09 next-gen-sequencing-annotations-ATACseq-KA.txt"
        self.df = pd.read_csv(self.input_metadata_file, sep="\t", index_col = "sample_id")
        logger.debug('self.df.shape{}'.format(self.df.shape))
        logger.debug('self.df.head{}'.format(self.df.head()))
# encoding="ISO-8859-1"

    def tearDown(self):
        print("Teardown")
        #delete the changed metadata files
        pass


    #def test_main_functional(self):
        """ usually when testing main we have to do that as a functional test - we pick a reasonable set of parameters for a 
        happy path and run through those to make sure that works
        """
    #    pm.main(None)


    #file exists
    #same file name

    def test_create_full_path(self):
        #checks if file exists
        metadata_dir = pm.create_full_path(self.metadata_subdir, self.metadata_file)
        logger.debug("Metadata file exists: {}".format(metadata_dir))
        self.assertTrue(os.path.exists(metadata_dir))
        #checks that the name of the path is correct
        self.assertEqual(metadata_dir, self.input_metadata_file)

    def test_load_metadata(self):
        # check the two dataframes are equal
        orig_metadata = pm.load_metadata(self.input_metadata_file)
        logger.debug("orig_metadata.shape: {}".format(orig_metadata.shape))
        self.assertEqual(orig_metadata.shape, self.df.shape)
        # print(self.df.shape)
        # print(self.df.size)
        # # print(self.df.columns)
        # print(self.df.iloc[6,4])
         

    def test_QC_sequence_removal(self):
        print("test_QC_sequence_removal")
        #samples_to_remove_from_metadata = [] # ["G8M54"]
        #edge case 1 --> remove 1 sample
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "c":["happy","sad", "laugh"], "d":["python","java", "html"]})
        logger.debug('\ntest_df:\n{}'.format(test_df))
        samples_to_remove = "c"
        updated_df = pm.QC_sequence_removal(test_df, samples_to_remove)
        logger.debug('\nupdated_df:\n{}'.format(updated_df))
        expected_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "d":["python","java", "html"]})
        logger.debug('\nexpected_df:\n{}'.format(expected_df))
        isEqual = updated_df == expected_df
        logger.debug('\nisEqual:\n{}'.format(isEqual))
        final_bool = isEqual.all()
        logger.debug('\nfinal_bool:\n{}'.format(final_bool))
        self.assertTrue(final_bool.all())
        #edge case 2 --> remove multiple samples
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "c":["happy","sad", "laugh"], "d":["python","java", "html"]})
        logger.debug('\ntest_df:\n{}'.format(test_df))
        samples_to_remove = ["a", "c"]
        updated_df = pm.QC_sequence_removal(test_df, samples_to_remove)
        logger.debug('\nupdated_df:\n{}'.format(updated_df))
        expected_df = pd.DataFrame({"b":["yellow","cherry", "book"], "d":["python","java", "html"]})
        logger.debug('\nexpected_df:\n{}'.format(expected_df))
        isEqual = updated_df == expected_df
        logger.debug('\nisEqual:\n{}'.format(isEqual))
        final_bool = isEqual.all()
        logger.debug('\nfinal_bool:\n{}'.format(final_bool))
        self.assertTrue(final_bool.all())

        #self.assertEqual(dataframe, )
        #set of samples to remove --> sample of index (check the intersection of the sets = empty)
        #dataframe pandas.index


    def test_check_number_group(self):
        #edge case 1 --> 3 replicates
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball", "apple", "crayon", "dog", "apple", "crayon", "dog"], "b":["yellow","cherry", "book", "yellow","cherry", "book", "yellow","cherry", "book"], "c":["happy","sad", "laugh", "apple", "crayon", "dog", "apple", "crayon", "dog"]})
        selected_cols = "b"
        all3s = pm.check_number_group(test_df, selected_cols, 3)
        logger.debug('\nall3s:\n{}'.format(all3s))
        checker = all3s == "3"
        logger.debug('\nchecker:\n{}'.format(checker))
        isequal = checker.all()
        logger.debug('\nisequal:\n{}'.format(isequal))
        self.assertTrue(isequal.all())
   
        #expected_replicates = pd.DataFrame()
        with self.assertRaises(pm.FhtbioinfpyPrepMetadataCheckValueReplicatesForGroups) as context:
            test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "c":["happy","sad", "laugh"]})
            selected_cols = "b"
            pm.check_number_group(test_df, selected_cols, 3, lambda : True)

        #edge case 3 --> replicates are not what were expected but user ignores exception error and continues
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["sample1","sample2", "sample3"], "c":["happy","sad", "laugh"]})
        selected_cols = "b"
        try:
            pm.check_number_group(test_df, selected_cols, 3, lambda : False)
        except pm.FhtbioinfpyPrepMetadataCheckValueReplicatesForGroups:
            self.fail("check_number_group() unexpectedly raised FhtbioinfpyPrepMetadataCheckValueReplicatesForGroups exception!")
        


    def test_add_string_columns(self):
        print("test_add_string_columns")
        test_df = pd.DataFrame({"pert_dose":[3, 4, 5], "pert_dose_unit": ["nM", "nM", "nM"], "pert_time":[10, 11, 12], "pert_time_unit": ["h", "h", "h"]})
        logger.debug('\ntest_df:\n{}'.format(test_df))
        new_df = pm.add_string_columns(test_df)
        logger.debug('\nnew_df:\n{}'.format(new_df))
        subset = False
        if set(['pert_dose_str','pert_time_str']).issubset(new_df.columns):
            subset = True
            #subset is true
        logger.debug('\nsubset:\n{}'.format(subset))
        #
        self.assertTrue(subset)
        self.assertEqual(new_df["pert_dose_str"][0],"3nM")
            
        #dataframe['pert_dose_str'] = dataframe['pert_dose'].astype(str) + dataframe['pert_dose_unit']
        #dataframe['pert_time_str'] = dataframe['pert_time'].astype(str) + dataframe['pert_time_unit']


    def test_verify(self):
        print("test_verify")

        #edge case 1 --> only one column in df is selected
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "c":["happy","sad", "laugh"], "d":["python","java", "html"]})
        col = "c"
        found = pm.verify(test_df, col)
        logger.debug('\nfound:\n{}'.format(found))
        self.assertTrue(found)

        #edge case 2 --> more than one column in df is selected
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "c":["happy","sad", "laugh"], "d":["python","java", "html"]})
        cols = ["a", "d"]
        found = pm.verify(test_df, cols)
        logger.debug('\nfound:\n{}'.format(found))
        self.assertTrue(found)

        with self.assertRaises(pm.FhtbioinfpyPrepMetadataVerifyColumnsForGroups) as context:
            #edge case 3 --> one column in df is selected but not found in df (print out exception)
            test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "c":["happy","sad", "laugh"], "d":["python","java", "html"]})
            cols = ["a", "f"] #f is a column not found in df
            pm.verify(test_df, cols)
        with self.assertRaises(pm.FhtbioinfpyPrepMetadataVerifyColumnsForGroups) as context:
            #edge case 4 --> more than one column not found in df (print out exception)
            test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "c":["happy","sad", "laugh"], "d":["python","java", "html"]})
            cols = ["a", "f", "h"] #columns f and h are not found in df
            pm.verify(test_df, cols)


        # break if it's a problem
        # check they're all in data frame column
        # set of data frame columns
        # set1 < set2

        # throw an exception --> add which columns aren't in the set
        

    def test_clean_selected_series(self):
        print("test_clean_columns")
        #edge case1 --> no cleaning necessary, dataframe before and after should be same
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"]})
        #logger.debug('\ntest_df:\n{}'.format(test_df))
        new_series = pm.clean_selected_series(test_df["a"])
        #logger.debug('\nnew_df:\n{}'.format(new_df))
        checker = new_series.equals(test_df["a"])
        #logger.debug('\nchecker:\n{}'.format(checker))
        self.assertTrue(checker)

        #edge case2 -->cleaning necessary
        test_series = pd.Series(["he-llo", "app)le", "basketball game"])
        logger.debug('\ntest_series:\n{}'.format(test_series))
        #expected_series = pd.Series({"a":["hello", "apple","basketball_game"]})
        #logger.debug('\nexpected_series:\n{}'.format(expected_series))
        new_series = pm.clean_selected_series(test_series)
        #new_series = cleaned["a"].squeeze()
        #new_series = new_series.reset_index(drop=True, inplace=True)
        logger.debug('\nnew_series:\n{}'.format(new_series))
        expected_series = pd.Series(["hello", "apple","basketball_game"])
        logger.debug('\nexpected_series:\n{}'.format(expected_series))
        isequal = new_series == expected_series
        logger.debug('\nisequal:\n{}'.format(isequal))
        self.assertTrue(isequal.all())

    def test_concatenateGroupDefColumns(self):
        print("test_concatenateGroupDefColumns")

        #edge case 1 --> one selected column
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "c":["happy","sad", "laugh"], "d":["python","java", "html"]})
        logger.debug('\ntest_df:\n{}'.format(test_df))
        selected_cols = ["c"]
        new_groups_df = pm.concatenateGroupDefColumns(test_df, selected_cols)
        logger.debug('\nnew_groups_df:\n{}'.format(new_groups_df))
        expected_groups_df = pd.DataFrame({"group":["happy", "sad", "laugh"]})
        logger.debug('\nexpected_groups_df:\n{}'.format(expected_groups_df))
        isEqual = new_groups_df == expected_groups_df
        logger.debug('\nisEqual:\n{}'.format(isEqual))
        final_bool = isEqual.all()
        logger.debug('\nfinal_bool:\n{}'.format(final_bool))
        self.assertTrue(final_bool.all())

        #edge case 2 --> more than one selected column
        test_df = pd.DataFrame({"a":["hello", "apple", "basketball"], "b":["yellow","cherry", "book"], "c":["happy","sad", "laugh"], "d":["python","java", "html"]})
        logger.debug('\ntest_df:\n{}'.format(test_df))
        selected_cols = ["a","c"]
        new_groups_df = pm.concatenateGroupDefColumns(test_df, selected_cols)
        logger.debug('\nnew_groups_df:\n{}'.format(new_groups_df))
        expected_groups_df = pd.DataFrame({"group":["hello_happy", "apple_sad", "basketball_laugh"]})
        logger.debug('\nexpected_groups_df:\n{}'.format(expected_groups_df))
        isEqual = new_groups_df == expected_groups_df
        logger.debug('\nisEqual:\n{}'.format(isEqual))
        final_bool = isEqual.all()
        logger.debug('\nfinal_bool:\n{}'.format(final_bool))
        self.assertTrue(final_bool.all())

        #selected_cols = ["pert_id"]
        #dataframe = pm.concatenateGroupDefColumns(self.df, selected_cols)
        #self.assertEqual(dataframe,)

    def test_add_groups_to_df(self):
        print("test_add_groups_to_df")
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
        final_bool = isEqual.all()
        logger.debug('\nfinal_bool:\n{}'.format(final_bool))
        self.assertTrue(final_bool.all())

    def test_add_experiment_id(self):
        print("test_add_experiment_id")
        experiment_id = "NS-21.0138"
        new_df = pm.add_experiment_id(self.df, experiment_id)
        if {'experiment_id'}.issubset(new_df):
            pass


    def test_create_R_Groups(self):
        print("test_create_R_Groups")

        # print('edge case 1 --> all start with a number / all are True')
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
        #print('edge case 2 --> none start with a number / all are False')
        test_df = pd.DataFrame({"a":["AB", "BC", "CD", "DE"], "b":["EF", "FG", "GH", "HI"]})
        logger.debug('test_df: \n{}\n'.format(test_df))
        new_df = pm.create_R_Groups(test_df, "b")
        logger.debug('new_Df: \n{}\n'.format(new_df))
        check1 = new_df.b.str[0:2]
        #logger.debug('\ncheck1:\n{}'.format(check1))
        check2 = check1 == "x_"
        #logger.debug('\ncheck2:\n{}'.format(check2))
        self.assertFalse(check2.all())
        
        # edge case 3 --> mix of strings that start with numbers and letters
        #print('edge case 3 --> mix of strings that start with numbers and letters')
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
    

    def test_save_to_csv(self):
        print("test_save_to_csv")

        #opening a file or directory --> make sure it gets cleaned up
        with tempfile.TemporaryDirectory(prefix = "fhtbioinfpy_test_prep_metadata") as wkdir:
            output_filepath = os.path.join(wkdir, "prep_metadata")
            pm.save_to_csv(output_filepath, self.df)
            self.assertTrue(os.path.exists(output_filepath))

        #make sure temporary folder is empty
        #make sure file exists
        #save to temporary folder
        #temp_dir
        


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()