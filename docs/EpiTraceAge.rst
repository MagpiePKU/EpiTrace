EpiTraceAge
-----------

Main function for age estimation. 

Input: 

- mat: a count matrix for scATAC/ATAC data, usually built with ArchR, SnapATAC or Signac etc. Row = peaks and col = cells. Specifically, each row (peak) is overlapping with ClockDML.- normalization_method: normalization method used in EpiTraceAge computation, choose between 'randomized', 'census', and 'blank'.     - ncores: number of cores used in parallel processing- size: depreciated parameter- subsamplesize: subsample cell number, default to 2000.    - select_minimal_percentage: selecting peaks covered in > this fraction of cells to compute dispersion. Default to 0.05. Suggest not to change unless necessary.            - select_size_of_dispersion: numbers of peaks selected for computing the similarity matrix, ranked by dispersion. Default to 3000. Suggest not to change unless necessary.  


Output: a data frame with columns: EpiTrace = epitrace age, Accessible_Loci = total accessibility on clock regions, cells=used cells.

 

  

