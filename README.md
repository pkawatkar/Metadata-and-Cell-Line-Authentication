# Metadata-and-Cell-Line-Authentication
prep_metadata, depmapID_lookup, cell_line_correlation, cell_line_authentication - summer internship 2022

    - prep_metadata module consists of modifications applied to a metadata file, describing samples in next_generation sequencing experiments, to regularize/standardize the metadata so that it can be used consistently with other analyses.

    - depmapID_lookup module is added on to the prep_metadata module to generate additional cleaned metadata about the cell line of each sample, which it then uses to lookup the reference ID used by the Broad Dependency Map and CCLE (Cancer Cell Line Encyclopedia).

    - cell_line_correlation module uses a "query" expression file and a "reference" expression file to calculate the corerlation values between each sample in the query with each cell line in the reference. This  was done by creating a subset of the input query_expr and input ref_expr that only has the intersection of the genes in each, and then using an external library(fast_corr) to efficiently calculate the correlation values.

    - cell_line_authentication module uses the correlation values from the cell_line_correlation module to find the most highly correlated cell line for each sample. It then checks whether the best correlated cell line matches the cell line that is annotated for the sample; if it does not, it prints an exception error of which of the samples do not match, as these samples may have incorrect metadata annotation or may indicate that the sequencing experiment had some technical issue.

