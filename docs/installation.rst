Installation
------------

EpiTrace requires R 4.0 or later. The current recommended environment is R 4.1.3 .


Install EpiTrace from Github_ using::

    library(devtools)   
    devtools::install_github('MagpiePKU/EpiTrace',dependencies=TRUE)  

If you encountered any issue in installation, please let us know by flagging up an issue `here in Github<https://github.com/MagpiePKU/EpiTrace/issues/new>`_.). 


Dependencies
^^^^^^^^^^^^

- dplyr
- tidyr
- RColorBrewer
- ggplot2
- Seurat (>=4.0) 
- Signac (>= 1.5.0)
- GenomicRanges
- plyranges (>= 1.0)
- WGCNA (>= 1.7)
- stringr
- easyLift
- parallel
- sva
- ccaPP
- HiClimR
- nnls
- ggtree
- ape
- reticulate
- ggpubr

**Session info** from developer environment::

	R version 4.2.0 (2022-04-22)
	Platform: aarch64-apple-darwin20 (64-bit)
	Running under: macOS 14.2

	Matrix products: default
	BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
	LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib

	locale:
	[1] C/UTF-8/C/C/C/C
	
	attached base packages:
	[1] parallel  stats4    stats     graphics  grDevices utils     datasets
	[8] methods   base

	other attached packages:
	 [1] ggpubr_0.4.0          reticulate_1.25       ape_5.7-1
	 [4] ggtree_3.4.0          nnls_1.4              HiClimR_2.2.1
	 [7] ccaPP_0.3.3           robustbase_0.95-0     pcaPP_2.0-1
	[10] sva_3.44.0            BiocParallel_1.30.2   genefilter_1.78.0
	[13] mgcv_1.8-40           nlme_3.1-157          easyLift_0.2.1
	[16] stringr_1.4.0         WGCNA_1.71            fastcluster_1.2.3
	[19] dynamicTreeCut_1.63-1 plyranges_1.16.0      GenomicRanges_1.48.0
	[22] GenomeInfoDb_1.34.9   IRanges_2.30.0        S4Vectors_0.34.0
	[25] BiocGenerics_0.42.0   ggplot2_3.3.6         RColorBrewer_1.1-3
	[28] tidyr_1.2.0           dplyr_1.1.4           Signac_1.6.0
	[31] sp_1.4-7              SeuratObject_4.1.0    Seurat_4.1.1
	[34] EpiTrace_0.0.0.9000


Again, if you run into issues in installation or using, do not hesitate to approach us or raise a `GitHub issue`_.

.. _Github: https://github.com/MagpiePKU/EpiTrace
.. _`Github issue`: https://github.com/MagpiePKU/EpiTrace/issues/new
