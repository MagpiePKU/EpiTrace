# EpiTrace
 Computing cell age with bulk and single-cell ATAC-seq data   

 EpiTrace takes an approximation approach to infer single cell age from single cell ATAC data. It infers single cell age by measuring the total opened reference genomic loci. On these loci, heterogeneity of chromatin accessibility decreases as the cell ages.   

 EpiTrace firstly algorithmically determine a set of tool genomic loci, on which the total chromatin accessibility (reads) shows maximal correlation to the total opened reference genomic loci. Then, the total chromatin accessibility on this set of tool genomic loci are used as an intermediate tool variable to approximate cell age.  
 
 EpiTrace documentation is now on `readthedocs`. 
 
 For descriptions, function references, and tutorials, visit https://epitrace.readthedocs.io 

 Maintainer: Zhang Yi <zhangyi@cimrbj.ac.cn>      

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

### ⚠️ Breaking Changes from v0.0.1

**NOT BACKWARDS COMPATIBLE with v0.0.1.x**

**Major API changes**:
- **Seurat v4 → v5**: EpiTrace now requires Seurat v5 (>= 5.0.0)
  - `GetAssayData(slot='data')` → `GetAssayData(layer='data')`
  - UMAP computation now uses UWOT (R package) instead of reticulate + Python
  - Some Seurat v4 workflows may need adjustments

**Removed dependencies**:
- ~~`ggtree`~~ → Use `ape::plot.phylo` instead for phylogeny
- ~~`easyLift`~~ → Use Bioconductor `easylift` instead for liftover

**If you're using EpiTrace v0.0.1**, you may need to update your code:
```r
# Old (v0.0.1)
library(ggtree)
plot <- ggtree(tree) + geom_tiplab()

# New (v0.0.2)
library(ape)
PlotEpiTracePhylogeny(phylo_result)  # New helper function
```

### System Requirements

**R Version**: >= 4.3.0 (recommended: >= 4.4.0)

**Key Dependencies** (updated in v0.0.2.0):
- **Seurat (>= 5.0.0)** ⚠️ Major upgrade from v4
- **SeuratObject** (Seurat v5)
- **Signac (>= 1.5.0)**
- **easylift** (Bioconductor) - NEW: replaces easyLift for liftover
- **ape** - NEW: replaces ggtree for phylogeny visualization
- **WGCNA (>= 1.7)**
- **GenomicRanges**
- **ggplot2**
- **dplyr**, **tidyr**, **RColorBrewer**
- **Matrix**, **matrixStats**, **sparseMatrixStats**
- **plyranges**
- **nnls**

**Full dependency list**: See [DESCRIPTION](DESCRIPTION)

**Platform**: Tested on macOS (ARM64/x86_64) and Linux

### Session Info

For reproducibility, here's the complete session info for EpiTrace v0.0.2.0:

```r
R version 4.4.2 (2024-10-31)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

Attached packages:
- EpiTrace_0.0.2.0
- Seurat_5.4.0              ⚠️ UPGRADED from v4.x
- SeuratObject_5.3.0         ⚠️ NEW in v5.x
- Signac_1.16.0
- easylift_1.7.0             ⚠️ NEW (replaces easyLift)
- ape_5.8                    ⚠️ NEW (replaces ggtree)
- WGCNA_1.73
- GenomicRanges_1.58.0
- ggplot2_3.5.1
- dplyr_1.1.4
```

**Key changes from previous versions**:
| Package | v0.0.1 | v0.0.2 | Breaking Change? |
|---------|--------|--------|------------------|
| Seurat  | v4.x    | v5.x   | ⚠️ YES - Major API changes |
| ggtree  | Used    | Removed | ⚠️ YES - Use `ape` instead |
| easyLift| Used    | Removed | ⚠️ YES - Use `easylift` instead |
| ape     | Optional| Required| ⚠️ YES - Now required |

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

