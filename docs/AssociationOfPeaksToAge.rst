AssociationOfPeaksToAge
-----------------------

Input: 

- epitrace_object (required): an age-determined EpiTrace object.  
- peakSetName: the "assay name" in the object that you would like to calculate age-x-peaks correlation with. Originally set to "peaks". But you can choose to calculate other associations such as "all". 
- epitrace_age_name: the vector name in meta.data for regression. Currently defaults to "EpiTraceAge_all". Can set this to "EpiTraceAge_iterative", for example. 
- subset_of_cells: a named string vector (cell name) for calculation. Default to NULL (no subsetting)
- epitrace_age_vector: if provided, regression would be performed on this vector instead of the named vector inside of EpiTrace object. 
- parallel: default to FALSE. Can set to TRUE to use multicore support. 

Output: a dataframe consisting name of each peak (rows), unscaled correlation coefficient (age-x-peak), and scaled correlation coefficient (age-x-peak).

The function is basically a wrapper of applying WGCNA::cor on a large matrix against a vector, with multicore support. 




