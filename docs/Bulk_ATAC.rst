Bulk ATAC example
-----------------

In this tutorial, we will show you a brief guide how to run EpiTrace analysis on bulk ATAC-seq, with a demo dataset from ``Corces 2016``.   
These demo data are from FACS-sorted pure hematopoietic cell populations and sequenced by ATAC-seq.   
Similar code works with single cell dataset.  

An important note is that EpiTrace age is reversed between bulk and single cell datasets. For a detailed technical explanation, please see ``Xiao 2022``. At the moment, please remember that higher age corresponds to **lower** EpiTrace age in bulk ATAC datasets, and **higher** EpiTrace age in single cell ATAC datasets. 

**Note this demo is an updated version from original *2022* version by using iterative updating algorithm**. We would suggest users to use the iterative updating algorithm regardless of whether it is bulk or single-cell ATAC-seq data.  


First of all, the input data for EpiTrace are: 

- a count matrix where rows are genomic loci, and columns are cells, and elements are ATAC-seq read numbers. 

- a data frame or GRanges object corresponds to the genomic loci (we call it "peak sets")   


You also need to choose a reference clock-like loci for age estimation. 
Currently, the package provides a human ClockDML dataset as reference. 
Alternative datasets could be chosen by you.  


Step 1. Initialize environment 
''''''''''''''''''''''''''''''

Load EpiTrace as::

    library(EpiTrace)

There are many other libraries should be loaded for this tutorial::

    library(dplyr)
    library(tidyr)
    library(readr)
    library(GenomicRanges)
    library(reshape2)
    library(openxlsx)
    library(ggplot2)
    library(Matrix)
    library(Seurat)
    library(SeuratObject)
    library(ggtree)
    library(EnsDb.Hsapiens.v86)
    library(patchwork)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(org.Hs.eg.db)
Step 2. Load data 
'''''''''''''''''
We load the demo human FACS-sorted pure blood cell data using::

    data("clock_gr_list")
    data("hematopoiesis.2018")


The package provided ``clock_gr_list`` is a list of GRanges, each corresponds to a set of clock-like DML (clockDML) in human.  

There are three sets of human clockDML that provided in `clock_gr_list`:  

`Mitosis`: 616 clock-like DML compiled from papers of ``Yang et al 2016``, ``Teschendorff 2020``, and ``Youn and Wang 2018``.     

`Chronology`: 125625 clock-like DML/DMR computed from genome-wide bisulfite capture sequencing data and a compilation of 450k methyl-array data, removing the ones that overlaps with Mitosis dataset.     

`solo_WCGW`: 1669720 solitary WCGW motif in partial methylated domain, from ``Zhou et al 2018``. This is just for comparison purpose and in general not used in age estimation.   


Step 3. Initialize peak sets 
''''''''''''''''''''''''''''
Peak sets are the annotation of ATAC peaks in your input dataset. Either a dataframe with columns chr,start,end or a GRanges object is fine. Here we provide you with an example GRanges ``inputpeak`` from the data ``hematopoiesis.2018``. Initaliztion of peak sets is necessary for following steps::   

    initiated_peaks <- Init_Peakset(inputpeak) 
    
Step 4. Initialize read count matrix
''''''''''''''''''''''''''''''''''''
We initiate a count matrix for EpiTrace by::   

    initiated_matrix <- Init_Matrix(cellname = colnames(inputmatrix),peakname = initiated_peaks$peakId, matrix = inputmatrix)
    

Step 5. Calculate EpiTrace Age
''''''''''''''''''''''''''''''
A one-shot function to calculate EpiTrace age using standard (human) ClockDML as reference:: 

    epitrace_obj_age_estimated <- EpiTraceAge_Convergence(initiated_peaks,inputmatrix,ref_genome = 'hg19',non_standard_clock = F,parallel = F,iterative_time = 10,Z_cutoff = 3,qualnum = 1)
    

Step 6. Compute the phylogenetic tree based on EpiTrace
'''''''''''''''''''''''''''''''''''''''''''''''''''''''
We pass the prepared epitrace object to estimate per-cell-cluster phylogenetic tree using RunEpiTracePhylogeny command.  

This generates a list which contains the assay (id of assay), tree (phylogenetic tree of the clusters), tree_plot (a ggtree object).  

If you would like to use other clustering of cells (say, single cell) you have to change the Idents of the object. Currently we tend to suggest use well annotated single cell clusters (or in this demo case, known FACS-sorted cell types).:: 

    Idents(epitrace_obj_age_estimated) <- epitrace_obj_age_estimated$celltype 
    phylotree_res_myeloid <- RunEpiTracePhylogeny(subset(epitrace_obj_age_estimated,celltype %in% c('HSC','MPP','CMP','GMP','Monocyte')))    chronology_tree_myeloid <- phylotree_res_myeloid[['iterative']][[2]]    chronology_tree_myeloid <- ape::root(chronology_tree_myeloid,outgroup='HSC')    plot(chronology_tree_myeloid)

Step 7. Compute association between regulatory region accessibility and cell age.
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
We pass the age-estimated object to AssociationOfPeaksToAge function to compute association result.  

Note that association is most meaningful for a sub-population or a lineage.:: 

    associated_res_myeloid <- AssociationOfPeaksToAge(subset(epitrace_obj_age_estimated,celltype %in% c('HSC','MPP','CMP','GMP','Monocyte')),epitrace_age_name = "EpiTraceAge_iterative")
    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
    associated_res_myeloid <- separate(associated_res_myeloid,col='locus',into=c('chr','start','end'),remove=F,convert=T)
    associated_res_myeloid_gr <- makeGRangesFromDataFrame(associated_res_myeloid)
    findOverlaps(associated_res_myeloid_gr,plyranges::reduce_ranges(c(clock_gr_list[[1]],clock_gr_list[[2]])))@from %>% unique -> peaks_overlap_with_clock
    associated_res_myeloid_gr$locus_type <- 'Others'
    associated_res_myeloid_gr$locus_type[peaks_overlap_with_clock] <- 'Clock-like DML'
    annotatePeak(associated_res_myeloid_gr, tssRegion=c(-2000, 500), TxDb=txdb,addFlankGeneInfo=F, flankDistance=50000,annoDb = "org.Hs.eg.db") -> associated_res_myeloid_gr_anno 
    as.data.frame(associated_res_myeloid_gr_anno@anno) -> associated_res_myeloid_gr_anno_df
    cbind(associated_res_myeloid, associated_res_myeloid_gr_anno_df %>% dplyr::select(SYMBOL,distanceToTSS,annotation,locus_type) ) -> associated_res_myeloid
    associated_res_myeloid$promoter_scaled_cor <- NA
    associated_res_myeloid$promoter_scaled_cor[grepl('Promoter', associated_res_myeloid$annotation)] <- scale(associated_res_myeloid $correlation_of_EpiTraceAge[grepl('Promoter', associated_res_myeloid$annotation)])
    associated_res_myeloid <- arrange(associated_res_myeloid,scaled_correlation_of_EpiTraceAge)
    (associated_res_myeloid) %>% na.omit() %>% tail(5) 

This shows the promoter of a monocyte-specific transcription factor MZB1 increases over age, in this lineage.  



**Reference**:   

Zhou W, et al. DNA methylation loss in late-replicating domains is linked to mitotic cell division. Nat Genet. 2018 Apr;50(4):591-602. doi: 10.1038/s41588-018-0073-4

Yang Z, et al. Correlation of an epigenetic mitotic clock with cancer risk. Genome Biol. 2016 Oct 3;17(1):205. doi: 10.1186/s13059-016-1064-3

Teschendorff AE. A comparison of epigenetic mitotic-like clocks for cancer risk prediction. Genome Med. 2020 Jun 24;12(1):56. doi: 10.1186/s13073-020-00752-3 

Youn A, Wang S. The MiAge Calculator: a DNA methylation-based mitotic age calculator of human tissue types. Epigenetics. 2018;13(2):192-206. doi: 10.1080/15592294.2017.1389361 

Horvath S. DNA methylation age of human tissues and cell types. Genome Biol. 2013;14(10):R115. doi: 10.1186/gb-2013-14-10-r115 

Corces MR, et al. Lineage-specific and single-cell chromatin accessibility charts human hematopoiesis and leukemia evolution. Nat Genet. 2016 Oct;48(10):1193-203. doi: 10.1038/ng.3646

Xiao Y, et al. Tracking single cell evolution via clock-like chromatin accessibility. BioRxiv. 2022. doi: 10.1101/2022.05.12.491736


