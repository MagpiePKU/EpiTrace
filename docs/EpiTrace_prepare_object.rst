EpiTrace_prepare_object
-----------------------

Input: 

- peakSet: a prepared "peak set" GRanges object. 
- matrix: a prepared read count matrix. 
- celltype: optional annotation for cell types. Usually not required. 
- sep_string: used to split the peak names. Usually as c(":","_","-"). 
- clock_gr_list: GRange list for reference clock-like loci. A list. 
- non_standard_clock: set to T when you need to use non-human, non-ClockDML genomic loci as reference. 
- qualnum: minimal peaks that a cell/sample should have. 
- remove_peaks_number:  minimal number of cells that a peak should be positive on. 
- ref_genome: when it is set to "hg19" or "hg38", the program can automatically use the ClockDML from package. Otherwise, you should provide your own reference clock-like loci. 
- run_reduction: set to F if you would not like Signac dimensionality reduction to run. 
- fn.k.param: KNN cluster parameter for Signac. 
- lsi_dim: LSI dimensions to be used in Signac dimensionality reduction. 
- min.cutoff: minimal cutoff for Signac::FindTopFeatures. Used when you ask the program to run dimensionality reduction. 

Output: an "epitrace" Seurat object to be used in downstream analysis. 

**NOTE** Do not use this function unless you intend to do the analysis step-by-step manually. The iterative updating function does not need object preparation at initial. The function is more internal and not recommended for application use. For application, we suggest directly using `EpiTraceAge_Convergence`. 
 

  

