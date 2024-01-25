

RunEpiTraceAge
--------------

Input: 

- epitrace_object: an prepared EpiTrace object. 
- parallel: set to T if you want multicore support. Requires ``parallel`` package. 
- ncores: threads used in parallel computation.  
- subsamplesize: cell number (randomly sampled) to be used in parallel computation. 

Output: an epitrace_object with cell ages determined. This is a Seurat object. In the object, most important results are:

- EpiTraceAge_xxxx. EpiTrace result from a particular (xxxx) reference clock-like loci. By default it runs for individual reference clock-like loci in `clock_gr_list`. The function is more internal and not recommended for application use. For application, we suggest directly using `EpiTraceAge_Convergence`. 



