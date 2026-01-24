.. role:: small
.. role:: smaller

Release Notes
=============


20250124 - v0.0.2.0
''''''''''''''''''''
**Major Release - Bug Fixes and Improvements**

**Breaking Changes:**
- Removed dependency on ``ggtree`` (use ``ape`` instead)
- Migrated from ``easyLift`` (GitHub) to ``easylift`` (Bioconductor)
- Updated for Seurat v5 compatibility

**Bug Fixes:**
- Fixed parse error with ref_genome parameter (Issue #15)
- Fixed hardcoded hg19 in EpiTraceAge_Convergence (Issue #12)
- Added offline mode support for UCSC chromosome download (Issue #3)
- Fixed vector mismatch warning in age calculation (Issue #19)
- Fixed final_cells undefined bug when run_reduction=FALSE

**New Features:**
- Added comprehensive documentation on bulk vs single-cell age interpretation (Issue #13)
- Added documentation explaining clock reference options (Issue #17)

**Dependencies Updated:**
- Removed: ``ggtree``, ``easyLift``
- Added: ``easylift`` (Bioconductor)
- Updated: Seurat v5 API compatibility (``slot`` → ``layer``)

**Notes:**
- All 9 fixes tested with bulk ATAC demo data (120 cells)
- EpiTrace now works offline without UCSC access
- Compatible with Bioconductor easylift for genome liftover


20240127
''''''''
Updated readthedocs introduction. 

Removed unused dependencies for now.  


20231011
''''''''
Included G-quadruplex reference set:


G-quadruplex (G4) harboring sites (pGQS) are provided for human (hg19), mouse (mm10), zebrafish (drerio10), fly (dm6), nematode (ce11) and yeast (sc3) genomes.
Source of G4 regions: GSE110582 (mm10, dm6, ce11, sc3), GSE187007 (hg19), or G4Hunter prediction (drerio10).


G4Hunter source: https://github.com/LacroixLaurent/G4HunterPaperGit


20221008
''''''''
Added frequently used non-human clockDML-like genomic regions as reference set for EpiTrace inference in: Drosophila, Zebrafish, and Mouse. 

20220910
''''''''
Added mouse ClockDML sites (from Zhou et al 2022, PMID 35873672) for mouse single cell EpiTrace inference.

20220611
''''''''
Added .rds files similar to the clock_gr_list data (for external loading).


Added %ni% function (taken from ArchR).

20220510
''''''''
First initialize version
