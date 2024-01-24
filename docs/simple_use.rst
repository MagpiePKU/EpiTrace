

API
===

Import EpiTrace as::

   library(EpiTrace)


The typical workflow consists of subsequent calls of
initializing the peak sets (``Init_Peakset``), 
initializing the read-count matrix (``Init_Matrix``), 
calculating the cell age using iterative updates (``EpiTraceAge_Convergence``). 
After age estimation, age-peak association (for identifying all clock-like loci: ``AssociationOfPeaksToAge``) and 
age-based cell phylogenic tree (``RunEpiTracePhylogeny``) could be calculated.  

If iterative updating of the reference clock-like loci is not required, 
we also provide functions to prepare a Seurat object which could be 
used by EpiTrace (``EpiTrace_prepare_object``), and 
to compute single-round EpiTrace age estimation from it (``RunEpiTraceAge``). 





