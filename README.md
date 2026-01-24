# EpiTrace
 Computing cell age with bulk and single-cell ATAC-seq data   

 EpiTrace takes an approximation approach to infer single cell age from single cell ATAC data. It infers single cell age by measuring the total opened reference genomic loci. On these loci, heterogeneity of chromatin accessibility decreases as the cell ages.   

 EpiTrace firstly algorithmically determine a set of tool genomic loci, on which the total chromatin accessibility (reads) shows maximal correlation to the total opened reference genomic loci. Then, the total chromatin accessibility on this set of tool genomic loci are used as an intermediate tool variable to approximate cell age.  
 
 EpiTrace documentation is now on `readthedocs`. 
 
 For descriptions, function references, and tutorials, visit https://epitrace.readthedocs.io 

 Maintainer: Zhang Yi <c.sinensis@gmail.com>      

### Installation

```r
if(!require(pak)){
    install.packages("pak")
}
library(pak)
pak::pkg_install('MagpiePKU/EpiTrace')
```

##### Development build (use at your own risk!)

```r
if(!require(pak)){
    install.packages("pak")
}
library(pak)
pak::pkg_install('MagpiePKU/EpiTrace@dev')
```

### System Requirements

**R Version**: >= 4.3.0

**Key Dependencies**:
- Seurat (>= 4.0)
- SeuratObject
- Signac (>= 1.5.0)
- easylift (Bioconductor) - for genome liftover
- ape - for phylogeny visualization (replaces ggtree)
- WGCNA (>= 1.7)
- GenomicRanges
- ggplot2

**Full dependency list**: See [DESCRIPTION](DESCRIPTION)

**Platform**: Tested on macOS (ARM64/x86_64) and Linux

### Session Info

For reproducibility, here's a typical session info for EpiTrace v0.0.2.0:

```r
R version 4.4.2 (2024-10-31)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.4 LTS

Other attached packages:
- EpiTrace_0.0.2.0
- Seurat_5.4.0
- Signac_1.16.0
- easylift_1.7.0
- ape_5.8
- GenomicRanges_1.58.0
```

### Changelog / Recent Updates

#### Version 0.0.2 (dev_2 branch)
Major bug fixes and improvements:

1. **Issue #15**: Fixed parse error - changed `%in%` to `==` for scalar ref_genome comparison
2. **Issue #12**: Removed hardcoded hg19 in `EpiTraceAge_Convergence` - now properly propagates ref_genome parameter
3. **Issue #3**: Added offline mode support - package now works without UCSC internet access
4. **Issue #19**: Fixed vector mismatch warning in age calculation by aligning vectors with `intersect()`
5. **ggtree removal**: Replaced with `ape::plot.phylo` for reduced dependencies
6. **easyLift → easylift**: Migrated from GitHub to Bioconductor for maintained package
7. **Chain file support**: Integrated built-in chain files from Bioconductor easylift
8. **Seurat v5 compatibility**: Updated deprecated `slot` parameter to `layer` in `GetAssayData`
9. **run_reduction fix**: Fixed `final_cells` undefined bug when `run_reduction=FALSE`

**Dependencies Updated**:
- Removed: `ggtree` (use `ape` instead)
- Removed: `easyLift` (GitHub, unmaintained)
- Added: `easylift` (Bioconductor, maintained)
- Updated: Seurat v5 API compatibility

### Citation
Xiao, Y., Jin, W., Ju, L. et al. Tracking single-cell evolution using clock-like chromatin accessibility loci. Nat Biotechnol (2024). https://doi.org/10.1038/s41587-024-02241-z

