
EpiTraceAge_Convergence
-----------------------

Main workhorse function in the package. 

Input: 

- peakSet: a prepared "peak set" GRanges object. 
- matrix: a prepared read count matrix. 
- celltype: optional annotation for cell types. Usually not required. 
- sep_string: used to split the peak names. Usually as c(":","_","-"). 
- clock_gr_list: GRange list for reference clock-like loci. A list. 
- non_standard_clock: set to T when you need to use non-default clock-like loci as reference. **NOTE** the standard clock-like loci is only useful when you work on human data. 
- qualnum: minimal peaks that a cell/sample should have. 
- remove_peaks_number:  minimal number of cells that a peak should be positive on. 
- ref_genome: when it is set to "hg19" or "hg38", the program can automatically use the clock-like loci from package. Otherwise, you should provide your own reference clock-like loci. 
- run_reduction: set to F if you would not like Signac dimensionality reduction to run. 
- fn.k.param: KNN cluster parameter for Signac. 
- lsi_dim: LSI dimensions to be used in Signac dimensionality reduction. 
- min.cutoff: minimal cutoff for Signac::FindTopFeatures. Used when you ask the program to run dimensionality reduction. 
- Z_cutoff: cutoff during iteration to select putative new clock-like loci. Defined as: include all loci with a scaled correlation coefficient peak x age > Z_cutoff. Usually set to 2~3. 
- mean_error_limit: limit of mean sample age differences between iterations. If the difference is smaller than this, iteration stops early. 
- iterative_time: maximal iteration time. 
- normalization_method: choose between "randomized", "census", and "blank". Default to "randomized". 
- select_minimal_percentage: minimal cell fraction that showing coverage on the selected peak. Default to 0.05. 
- select_size_of_dispersion: number of peaks selected for computing the similarity matrix, ranked by dispersion. 


Output: a Seurat object with ages annotated in its metadata. 

Main informations added to the meta.data of the Seurat objects are: 

- EpiTraceAge_Clock_initial: the age estimated using only reference clock-like loci 
- EpiTraceAge_iterative: the age estimated using reference + newly-included (putative sample-specific) clock-like loci 


