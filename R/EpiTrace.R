# EpiTrace: package for computing cell age and phylogeny from ATAC-seq data
# Zhang Yi, <synapse@pku.edu.cn>


# library(glmnet)
# library(dplyr)
# library(tidyr)
# library(readr)
# library(dplyr)
# library(ggplot2)
# library(reshape2)
# library(stringr)
# library(ggbeeswarm)
# library(ggsci)
# library(ggpubr)
# library(openxlsx)
# library(ArchR)
# library(ggtree)


#' Init_Peakset: function for initiating the peak set object for EpiTrace.
#' @title Init_Peakset
#'
#' @description function for converting an peakset to useful object in EpiTrace.
#'
#' @details Init_Peakset(peakset)
#'
#' @param peakset either a GRanges object, or a data frame (at least) with columns named 'chr','start',and 'end'.
#' @return a GRanges object with sequenced peak identifier (peakId).
#' @export
#' @examples Init_Peakset(peakset)
#'


Init_Peakset <- function(peakSet){
  if(class(peakSet)[1] %in% 'GRanges'){
    peakSet <- peakSet
  }else{
    peakSet <- peakSet %>% makeGRangesFromDataFrame()
  }
  peakSet$peakId <- 1:length(peakSet)
  return(peakSet)
}

#' Init_Matrix: function for initiating the matrix for EpiTrace.
#' @title Init_Matrix
#'
#' @description function for converting and checking coherence of an input matrix in EpiTrace.
#'
#' @details Init_Matrix(cellname,peakname,matrix)
#'
#' @param peakname a vector
#' @param cecllname a vector
#' @param matrix the input matrix
#'
#' @return a matrix with fixed cell name/peak name
#' @export
#' @examples Init_Matrix(cellname,peakname,matrix)
#'



Init_Matrix <- function(cellname,peakname,matrix){
  if(length(cellname)!=ncol(matrix)){
    message('cell names number is not similar to cells in matrix ')
    stop()
  }
  if(length(peakname)!=nrow(matrix)){
    message('peak names number is not similar to peaks in matrix ')
    stop()
  }
  colnames(matrix) <- cellname
  rownames(matrix) <- peakname
  return(matrix)
}


#' EpiTrace_prepare_object: wrapper function for preparing input data matrix for EpiTrace.
#' @title EpiTrace_prepare_object
#'
#' @description wrapper function for preparing input data matrix for EpiTrace.
#'
#' @details EpiTrace_prepare_object (peakSet,matrix,celltype=NULL,min.cutoff=50,lsi_dim=2:50,fn.k.param=21,ref_genome='hg38',sep_string= c(":", "-"),clock_gr_list=clock_gr_list,non_standard_clock=F)
#' @details Note: SC/Bulk could not be run at the same time as the normalization/filtering process kills single cells in the face of bulk data
#'
#' @param matrix input count matrix of ATAC data, rows=GRanges and cols=samples/single cells.
#' @param peakSet a GenomicRanges object corresponding to the rows of count matrix.
#' @param min.cutoff minimal cutoff for Signac variableFeature calling
#' @param lsi_dim the dimensionalities used in LSI for Signac
#' @param fn.k.param the k parameter used in FindNeighbours
#' @param ref_genome hg38 or hg19. Currently only support these two.
#' @param sep_string the separation string for row names in input matrix to generate ranges. for example, 'chr1:1-2' is c(':','-')
#' @param clock_gr_list the clockDML set, is a list of reference GRanges, each corresponds to a set of clock-like DML or DMR
#' @param non_standard_clock whether is not using non-standard reference clockDML sets.
#' @param run_reduction whether run reduction or not. Could be time saving to not run it if you only need age derivation. 
#'
#' @return a seurat object with all input peaks, as well as assays named 'x' for each clock set.
#' @export
#' @examples


EpiTrace_prepare_object <- function(peakSet,matrix,celltype=NULL,min.cutoff=50,lsi_dim=2:50,fn.k.param=21,ref_genome='hg38',sep_string= c(":", "-"),clock_gr_list=clock_gr_list,non_standard_clock=F,qualnum=10,run_reduction=T){
  # test
  # initiated_peaks -> peakSet
  # singleCell_Hemapoietic_dat -> matrix
  # cellnames_singlecell$celltype -> celltype
  # ref_genome = 'hg38'
  # non_standard_clock = T
  # clock_gr_list = target_clock_region_for_sc
  # min.cutoff=50
  # lsi_dim=2:50
  # fn.k.param=21
  # sep_string= c(":", "-")
  
  
  # 0. try error catch
  if(!is.null(celltype)){
    if(length(celltype) != ncol(matrix)){
      message('Input cell type vector does not match the cell# in input data')
      stop()
    }
  }
  if(!is.null(peakSet)){
    if(length(peakSet) != nrow(matrix)){
      message('Input peak GRanges# does not match the ranges# in input data')
      stop()
    }
  }
  if(ref_genome != 'hg19' & ref_genome != 'hg38'){
    message('ref genome is neither hg19 nor hg38')
    # stop()
  }
  if(non_standard_clock %in% F){
    if(sum(names(clock_gr_list) %in% c('Mitosis','Chronology','solo_WCGW'))!=3){
      message('ref clock list is said to be standard but actually not')
      stop()
    }else{
      message('ref clock list is set to be standard (Homo sapiens, hg19)')
      standard_clock = T
    }
  }else{
    message('ref clock list is not standard. Please make sure the input data, peak set and clock set are in similar reference genome.')
    standard_clock = F
  }
  message('Input peakset is set to be ',ref_genome)
  
  # 0.5. Filter the input matrix to a better one (try removing zero expression cells, and zero expression peaks)
  which(colSums(matrix)==0) %>% names -> remove_cells
  which(rowSums(matrix)<10 | rowSums(matrix>0)<10) %>% names -> remove_peaks
  logical_cell_vec <- colnames(matrix) %in% remove_cells
  logical_peak_vec <- rownames(matrix) %in% remove_peaks
  #peakSet,matrix,celltype
  matrix -> original_matrix
  if(length(remove_cells)>0){
    if(!is.null(celltype)){
      celltype -> original_celltype
      celltype <- celltype[!logical_cell_vec]
    }
    matrix <- matrix[,!logical_cell_vec]
  }
  
  if(length(remove_peaks)>0){
    peakSet -> original_peakSet
    peakSet <- peakSet[!logical_peak_vec,]
    matrix <- matrix[!logical_peak_vec,]
  }
  if(!is.null(celltype)){
    names(celltype) <- colnames(matrix)
  }
  
  
  # 1. prepare Clock region intersection to input peak set.
  Overlap_Input_with_Clock(peakSet_generanges=peakSet,clock_gr_list=clock_gr_list,ref=ref_genome) -> overlap_result
  overlap_list_of_list <- overlap_result$overlap_list_of_list
  if(standard_clock & (ref_genome %in% 'hg38')){
    clock_gr_list[['Mitosis']] %>% easyLift::easyLiftOver('hg19_hg38') -> mitosis_gr
    clock_gr_list[['Chronology']] %>% easyLift::easyLiftOver('hg19_hg38') -> chronology_gr
    plyranges::bind_ranges(mitosis_gr,chronology_gr) %>% reduce()  -> target_clock_gr
    result_clock_gr_list <- list('MitosisClock'=mitosis_gr,'ChronologyClock'=chronology_gr,'AllClock'=target_clock_gr)
  }
  if(standard_clock & (ref_genome %in% 'hg19')){
    clock_gr_list[['Mitosis']] -> mitosis_gr
    clock_gr_list[['Chronology']] -> chronology_gr
    plyranges::bind_ranges(mitosis_gr,chronology_gr) %>% reduce()  -> target_clock_gr
    result_clock_gr_list <- list('MitosisClock'=mitosis_gr,'ChronologyClock'=chronology_gr,'AllClock'=target_clock_gr)
  }
  if(non_standard_clock){
    result_clock_gr_list <- clock_gr_list
  }
  
  
  # 2. prepare the chromatin assay with Signac
  Signac::CreateChromatinAssay(matrix, sep = sep_string,
                               genome = ref_genome,ranges=peakSet) -> chrom_assay
  # 3. prepare the Seurat object
  tempdf <- Seurat::CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks"
  )
  # 4. add clock sets
  # modified
  # add try-catch to remove zero-expression cells. leaving only those cells OK for all clock measurements. Of course this could be modified in the future but at the moment let it be.
  for (x in names(result_clock_gr_list)){
    message('add ',x)
    matrix[findOverlaps(peakSet,result_clock_gr_list[[x]])@from %>% unique,] -> mtx2
    message('good quality cells ',sum(colSums(mtx2,na.rm=T)>=qualnum),' and peaks ',sum(rowSums(mtx2,na.rm=T)>=qualnum))
    dim(mtx2) -> dimcurr
    p=0
    count = 0
    while(p==0){
      mtx2 <- mtx2[rowSums(mtx2,na.rm=T)>=qualnum,colSums(mtx2,na.rm=T)>=qualnum]
      dim(mtx2) -> dimnext
      if((dimcurr[1]-dimnext[1])==0){
        p=1
      }
      count = count + 1
      # print(dimcurr[1]-dimnext[1]) # remove this debug tag
    }
    tempdf$cell <- rownames(tempdf@meta.data)
    tempdf <- subset(tempdf,cell %in% colnames(mtx2))
    tempdf[[x]] <- Seurat::CreateAssayObject(counts = mtx2[,tempdf$cell],min.cells = 0,min.features = 0,check.matrix = F)
  }
  # 5. perform dimensionality reduction, feature selection, etc, on this final object
  if(run_reduction==T){
    tempdf <- Signac::RunTFIDF(tempdf)
    tempdf <- Signac::FindTopFeatures(tempdf, min.cutoff = min.cutoff)
    tempdf <- Signac::RunSVD(tempdf)
    tempdf <- Seurat::RunUMAP(object = tempdf, reduction = 'lsi', dims = lsi_dim)
    tempdf <- Seurat::FindNeighbors(object = tempdf, reduction = 'lsi', dims = lsi_dim,k.param = fn.k.param)
    tempdf <- Seurat::FindClusters(object = tempdf, verbose = FALSE, algorithm = 3)
    
    # define identity
    tempdf@meta.data %>% rownames -> final_cells
    if(!is.null(celltype)){
      tempdf@meta.data$celltype <- celltype[final_cells]
      Idents(tempdf) <- celltype[final_cells]
    }else{
      tempdf@meta.data$celltype <- tempdf@meta.data$seurat_clusters
      Idents(tempdf) <- tempdf$seurat_clusters
    }
  }else{
    if(!is.null(celltype)){
      tempdf@meta.data$celltype <- celltype[final_cells]
      Idents(tempdf) <- celltype[final_cells]
    }
  }

  return(tempdf)
}




#' RunEpiTraceAge: wrapper function for computing EpiTrace age for input Seurat object with scATAC/bulkATAC data.
#' @title RunEpiTraceAge
#'
#' @description wrapper function for computing EpiTrace age for input scATAC/bulkATAC data.
#'
#' @details RunEpiTraceAge(epitrace_object)
#' @details Note: SC/Bulk could not be run at the same time as the normalization/filtering process kills single cells in the face of bulk data
#'
#' @param epitrace_object a seurat object built by EpiTrace_prepare_object
#' @return a Seurat_Object with added meta.data columns EpiTraceAge_Mitotic, EpiTraceAge_Chronological, EpiTraceAge_all, Accessibility_Mitotic, Accessibility_Chronological, Accessibility_all. Where the 'Mitotic' corresponds to EpiTraceAge computed with mitosis-associated DML and 'Chronological' corresponds to non-mitosis, clock-like DML. 'all' is the EpiTraceAge computed by both sets.
#' @export
#' @examples
#'


RunEpiTraceAge <- function(epitrace_object,parallel=F,ncores=20,subsamplesize=2000){
  # test
  # epitrace_object <- epitrace_obj
  availableAssays <- SeuratObject::Assays(epitrace_object)
  availableAssays_non_peak <- availableAssays[availableAssays!='peaks']
  epitrace_object$cell <- rownames(epitrace_object@meta.data)
  mtx_list <- lapply(availableAssays_non_peak,function(x){
    DefaultAssay(epitrace_object) <- x
    Seurat::GetAssayData(epitrace_object,slot='data')
  })
  names(mtx_list) <- availableAssays_non_peak
  lapply(names(mtx_list),function(x){
    mtx_list[[x]] -> testmtx
    if(sum(Matrix::rowSums(testmtx)>0,na.rm=T)==0 & sum(Matrix::colSums(testmtx)>0,na.rm=T)==0){
      message(paste('no data overlapping for clock reference set',x))
      stop()
    }else{
      message('Estimating Age by EpiTraceAge...')
      if(parallel==F){
        EpiTraceAge(mat = testmtx[Matrix::rowSums(testmtx)>0,Matrix::colSums(testmtx)>0])
      }else{
        EpiTraceAge(mat = testmtx[Matrix::rowSums(testmtx)>0,Matrix::colSums(testmtx)>0],ncores = ncores,parallel=T,subsamplesize=subsamplesize)
      }
    }
  }) -> res_list
  names(res_list) <- names(mtx_list)
  count = 0
  for (x in names(res_list)){
    count = count + 1
    res_list[[x]] -> df1
    colnames(df1)[2:4] <- paste0(c('EpiTraceAge_','Accessibility_','AccessibilitySmooth_'),x)
    if(count == 1){
      returndf = df1
    }else{
      returndf <- left_join(returndf,df1)
    }
  }
  epitrace_object@meta.data$cell <- rownames(epitrace_object@meta.data)
  left_join(epitrace_object@meta.data %>% as.data.frame(),returndf) -> cell_res
  rownames(cell_res) <- rownames(epitrace_object@meta.data)
  epitrace_object@meta.data <- cell_res
  return(epitrace_object)
}




#' RunEpiTracePhylogeny: wrapper function for computing EpiTrace phylogeny for input Seurat object with scATAC/bulkATAC data.
#' @title RunEpiTracePhylogeny
#'
#' @description wrapper function for computing EpiTrace phylogeny for input scATAC/bulkATAC data.
#'
#' @details RunEpiTracePhylogeny(epitrace_object)
#'
#' @param epitrace_object a seurat object prepared by EpiTrace_prepare_object
#' @param min.cutoff min.cutoff used in Signac::FindTopFeatures analysis, during re-selecting the clock variable features
#' @return a list, in which each element corresponds to a clock dataset in the original object. Each element is a list: clock = the name of clockDML set, tree = phylogenetic tree for 'clusters of cells' defined by the given 'idents' of the input seurat object, and tree_plot = a ggtree plot of the tree.
#' @export
#' @examples
#'


RunEpiTracePhylogeny <- function(epitrace_object,min.cutoff=50,run_reduction=T){
  availableAssays <- SeuratObject::Assays(epitrace_object)
  epitrace_object$cell <- rownames(epitrace_object@meta.data)
  color_celltype <- grDevices::colorRampPalette(colors = RColorBrewer::brewer.pal(9,"Set1"))(length(unique(Idents(epitrace_object))))
  names(color_celltype) <- unique(Idents(epitrace_object))
  lapply(availableAssays,function(assayid){
    tryCatch({
      if(DefaultAssay(epitrace_object) != assayid){
        DefaultAssay(epitrace_object) <- assayid
        if(run_reduction==T){
          epitrace_object <- Signac::RunTFIDF(epitrace_object)
          epitrace_object <- Signac::FindTopFeatures(epitrace_object, min.cutoff = min.cutoff)
          epitrace_object <- Signac::RunSVD(epitrace_object)
        }
      }
      obj_clock <- BuildClusterTree(object = epitrace_object,verbose = T,assay = assayid)
      data.tree_clock <- Tool(object = obj_clock, slot = "BuildClusterTree")
      ggtree::ggtree(data.tree_clock,layout='rectangular',ladderize = FALSE)  + geom_tiplab(aes(color=label),size=5,offset=10) + scale_color_manual(values=color_celltype)  -> tree_plot_clock
      xmax <- (tree_plot_clock$data$`branch.length` %>% max(na.rm=T)) * 1.4
      tree_plot_clock <- tree_plot_clock + xlim(c(NA,xmax))
      result <- list(assay=assayid,tree=data.tree_clock,tree_plot=tree_plot_clock)
      return(result)
    },error=function(e){message('failed for ',assayid)})
  }) -> returnlist
  names(returnlist) <- availableAssays
  return(returnlist)
}



#' Overlap_Input_with_Clock: function for overlapping the input peak set object to clockDML for EpiTrace.
#' @title Overlap_Input_with_Clock
#'
#' @description function for converting an input peakset to a list of peaks that overlapping with ClockDML in EpiTrace.
#'
#' @details Overlap_Input_with_Clock(peakSet_generanges,clock_gr_list=clock_gr_list,ref=ref)
#'
#' @param peakSet_generanges the input peak set, a GRanges object
#' @param clock_gr_list a list of GRanges, each corresponds to one Clock-DML set. Standard reference is made accompanying the package with 'Mitosis','Chronological','solo_WCGW' and 'all' (Mitosis+Chronological). These are all in hg19.
#' @ref The reference genome of input GRanges. 'hg38' or 'hg19'.
#'
#' @return a list with: overlap_list_of_list corresponds to Ids of peaks overlapping with ClockDML, and a dataframe overlap_df which tells the stats of overlapping. Currently the function would NOT tell a problem when there is no overlapping peak is found.
#' @export
#' @examples Overlap_Input_with_Clock(peakSet_generanges=test_Peakset,clock_gr_list=clock_gr_list,ref='hg38')
#'

Overlap_Input_with_Clock <- function(peakSet_generanges,clock_gr_list=clock_gr_list,ref=ref){
  if(ref %in% 'hg19'){
    # do nothing
  }
  if(ref %in% 'hg38'){
    lapply(names(clock_gr_list),function(x){
      easyLift::easyLiftOver(clock_gr_list[[x]],'hg19_hg38') -> temp
      return(temp)
    }) -> clock_gr_list_new
    names(clock_gr_list_new) <- names(clock_gr_list)
    clock_gr_list_new -> clock_gr_list
  }
  ref_gr_list <- list('Input'=peakSet_generanges)
  lapply(names(clock_gr_list),function(x){
    lapply(names(ref_gr_list),function(y){
      findOverlaps(ref_gr_list[[y]],clock_gr_list[[x]])
    }) -> returnlist
    names(returnlist) <- paste0(names(ref_gr_list),'_vs_',x)
    returnlist
  }) -> overlap_list_of_list
  names(overlap_list_of_list) <- names(clock_gr_list)
  lapply(overlap_list_of_list,function(x) lapply(x,function(y) y@to %>% unique %>% length()) %>% unlist) %>% unlist -> overlap_num
  overlap_num <- data.frame(name=names(overlap_num),overlap_num=overlap_num)
  overlap_num <- separate(overlap_num,col=1,into=c('Clock_panel','Comparison'),remove=F,sep='\\.')
  overlap_num$target <- gsub('_vs.+','',overlap_num$Comparison)
  data.frame(Clock_panel=names(clock_gr_list),Size=lapply(clock_gr_list,function(x) nrow(as.data.frame(x)))%>%unlist)-> clock_size
  overlap_df <- left_join(overlap_num,clock_size)
  overlap_df$overlap_fraction <- overlap_df$overlap_num/overlap_df$Size
  return(list(overlap_list_of_list=overlap_list_of_list,overlap_df=overlap_df))
}


#' EpiTraceAge: function for computing cell age with EpiTrace for input scATAC/bulkATAC data.
#' @title EpiTraceAge
#'
#' @description function for computing cell age with EpiTrace for input scATAC/bulkATAC data.
#'
#' @details EpiTraceAge(mat)
#' @details Note: SC/Bulk could not be run at the same time as the normalization/filtering process kills single cells in the face of bulk data
#' @details This is adapted from CytoTrace method.
#'
#' @param mat a count matrix for scATAC/ATAC data, usually built with ArchR, SnapATAC or Signac etc. Row = peaks and col = cells. Specifically, each row (peak) is overlapping with ClockDML.#'
#'
#' @return a data frame with EpiTrace = epitrace age, Accessible_Loci = total accessibility on clock regions,cells=used cells).
#' @export
#' @examples


EpiTraceAge_backup <- function(mat){
  require(Matrix)
  bad_peaks <- rowSums(mat) == 0
  num_bad_peaks <- length(which(bad_peaks == TRUE))
  mat <- mat[!bad_peaks,]
  mat_log <- log(mat+1,2)
  mat_log <- data.matrix(mat_log)
  counts <- rowSums(mat>0)
  reads <- rowSums(mat)
  normalized_mat <- t(t(mat)*counts/reads)
  normalized_mat <- log(normalized_mat+1,2)
  accessible_loci <- rowSums(normalized_mat > 0)
  per_cell_accessible_loci <- colSums(normalized_mat > 0 )
  peaks_expressed_in_n_cell <- rowSums(normalized_mat>0)
  peaks_expressed_in_n_cell>0.05*ncol(mat) -> good_peaks_for_selecting_loci
  mtx_with_better_peaks <- normalized_mat[good_peaks_for_selecting_loci,]
  if(class(mtx_with_better_peaks)[[1]]=='dgeMatrix'){
    dispersion <- matrixStats::rowVars(mtx_with_better_peaks %>% as.matrix())/rowMeans(mtx_with_better_peaks)
  }else{
    dispersion <- sparseMatrixStats::rowVars(mtx_with_better_peaks)/rowMeans(mtx_with_better_peaks)
  }
  size_variable_loci <- 3000
  cutoff <- sort(dispersion) %>% tail(min(length(dispersion),size_variable_loci)) %>% head(1) %>% as.numeric()
  normalized_variable_mat <- normalized_mat[dispersion >= cutoff,]
  normalized_variable_mat <- normalized_variable_mat[,colSums(normalized_variable_mat)>0]
  WGCNA::cor(normalized_variable_mat) -> cor_normalized_variable_mat
  mean(as.vector(cor_normalized_variable_mat)) -> mean_number
  diag(cor_normalized_variable_mat) <- 0
  cor_normalized_variable_mat[cor_normalized_variable_mat<0] <- 0
  cor_normalized_variable_mat[cor_normalized_variable_mat<mean_number] <- 0
  censored_cor_normalized_variable_mat <- cor_normalized_variable_mat/rowSums(cor_normalized_variable_mat)
  censored_cor_normalized_variable_mat[rowSums(cor_normalized_variable_mat)==0,] <- 0
  final_cells <- colnames(censored_cor_normalized_variable_mat)
  mat2 <- normalized_mat[,final_cells]
  per_cell_accessible_loci <- per_cell_accessible_loci[final_cells]
  correlation_of_each_peak_to_per_cell_accessible_loci_counts <- WGCNA::cor(x=t(mat2),y=per_cell_accessible_loci,drop = T)
  names(sort(correlation_of_each_peak_to_per_cell_accessible_loci_counts)%>%tail(min(length(correlation_of_each_peak_to_per_cell_accessible_loci_counts),200))) -> best_loci
  accessible_loci_counts_smoothened <- colMeans(mat2[best_loci,])
  nnls_res <- nnls::nnls(censored_cor_normalized_variable_mat,accessible_loci_counts_smoothened)
  nx1_mat <- censored_cor_normalized_variable_mat %*% nnls_res$x
  stepcount = 1
  delta = 1
  extension_parameter=0.99
  original_nx1_mat <- nx1_mat
  current_nx1_mat <- nx1_mat
  next_nx1_mat <- nx1_mat
  while(stepcount < 50000 & delta >= 1e-9){
    stepcount = stepcount + 1
    next_nx1_mat <- extension_parameter*(censored_cor_normalized_variable_mat%*% current_nx1_mat) + (1-extension_parameter) * original_nx1_mat
    delta <- mean(abs(next_nx1_mat-current_nx1_mat))
    current_nx1_mat <- next_nx1_mat
  }
  epitrace_ranked <- rank(current_nx1_mat)
  epitrace_norm <- (epitrace_ranked-min(epitrace_ranked))/(max(epitrace_ranked)-min(epitrace_ranked))
  data.frame(
    cell=final_cells,
    EpiTrace = epitrace_norm,
    Accessible_Loci = per_cell_accessible_loci,
    Accessible_Loci_Smooth = accessible_loci_counts_smoothened
  ) -> returndf
  return(returndf)
}

# subsample, switch between parallel/nonparallel 
							  
EpiTraceAge <- function(mat,ncores=20,parallel=F,size=2000,batch=NULL,subsamplesize=2000){
  #test
  # mat <- mat_test
  # ncores=20
  # parallel=F
  # size=2000
  # batch=NULL
  # subsamplesize=2000
  
  #Subsample routine
  if(parallel == FALSE){
    size <- ncol(mat)
  } else {
    size <- min(ncol(mat),subsamplesize)
  }
  chunk <- round(ncol(mat)/size)
  subsamples <- split(1:ncol(mat), sample(factor(1:ncol(mat) %% chunk)))
  require(Matrix)
  bad_peaks <- rowSums(mat) == 0
  num_bad_peaks <- length(which(bad_peaks == TRUE))
  mat <- mat[!bad_peaks,]
  chunk <- round(ncol(mat)/size)
  subsamples <- split(1:ncol(mat), sample(factor(1:ncol(mat) %% chunk)))
  mat_origin <- mat
  batch_origin <- batch
  sub_batch_function <- function(subsample){
    #Checkpoint: Sequencing depth normalization
    mat <- mat_origin[,subsample]
    batch <- batch_origin[subsample]
    mat_log <- log(mat+1,2)
    mat_log <- data.matrix(mat_log)
    counts <- rowSums(mat>0,na.rm=T)
    reads <- rowSums(mat,na.rm=T)
    normalized_mat <- t(t(mat)*counts/reads)
    normalized_mat <- log(normalized_mat+1,2)
    normalized_mat[is.na(normalized_mat)] <- 0
    accessible_loci <- rowSums(normalized_mat > 0)
    per_cell_accessible_loci <- colSums(normalized_mat > 0)
    peaks_expressed_in_n_cell <- rowSums(normalized_mat>0)
    peaks_expressed_in_n_cell>0.05*ncol(mat) -> good_peaks_for_selecting_loci
    mtx_with_better_peaks <- normalized_mat[good_peaks_for_selecting_loci,]
    if(class(mtx_with_better_peaks)[[1]]=='dgeMatrix'){
      dispersion <- matrixStats::rowVars(mtx_with_better_peaks %>% as.matrix())/rowMeans(mtx_with_better_peaks)
    }else{
      dispersion <- sparseMatrixStats::rowVars(mtx_with_better_peaks)/rowMeans(mtx_with_better_peaks)
    }
    size_variable_loci <- 3000
    cutoff <- sort(dispersion) %>% tail(min(length(dispersion),size_variable_loci)) %>% head(1) %>% as.numeric()
    normalized_variable_mat <- normalized_mat[names(dispersion)[dispersion >= cutoff],]
    normalized_variable_mat <- normalized_variable_mat[,colSums(normalized_variable_mat)>0]
    WGCNA::cor(normalized_variable_mat) -> cor_normalized_variable_mat
    mean(as.vector(cor_normalized_variable_mat)) -> mean_number
    diag(cor_normalized_variable_mat) <- 0
    cor_normalized_variable_mat[cor_normalized_variable_mat<0] <- 0
    cor_normalized_variable_mat[cor_normalized_variable_mat<mean_number] <- 0
    censored_cor_normalized_variable_mat <- cor_normalized_variable_mat/rowSums(cor_normalized_variable_mat)
    censored_cor_normalized_variable_mat[rowSums(cor_normalized_variable_mat)==0,] <- 0
    final_cells <- colnames(censored_cor_normalized_variable_mat)
    mat2 <- normalized_mat[,final_cells]
    per_cell_accessible_loci <- per_cell_accessible_loci[final_cells]
    list(mat2 = mat2,counts = counts, censored_cor_normalized_variable_mat = censored_cor_normalized_variable_mat,per_cell_accessible_loci=per_cell_accessible_loci) -> returnlist
    return(returnlist)
  }
  if(parallel==T){
    batches <- parallel::mclapply(subsamples, mc.cores = min(chunk, ncores), function(subsample_x){
      sub_batch_function(subsample=subsample_x)
    })
  }else{
    batches <- lapply(subsamples, function(subsample_x){
      sub_batch_function(subsample=subsample_x)
    })
  }
  mat2 <- do.call(cbind, lapply(batches, function(x) x$mat2))
  counts <- do.call(c, lapply(batches, function(x) x$counts))
  per_cell_accessible_loci <- lapply(batches,function(x) x$per_cell_accessible_loci) %>% unlist()
  correlation_of_each_peak_to_per_cell_accessible_loci_counts <- WGCNA::cor(x=t(mat2),y=per_cell_accessible_loci,drop = T)
  correlation_of_each_peak_to_per_cell_accessible_loci_counts[is.na(correlation_of_each_peak_to_per_cell_accessible_loci_counts)] <- 0
  names(sort(correlation_of_each_peak_to_per_cell_accessible_loci_counts)%>%tail(min(length(correlation_of_each_peak_to_per_cell_accessible_loci_counts),200))) -> best_loci
  accessible_loci_counts_smoothened <- colMeans(mat2[best_loci,])
  samplesize <- unlist(lapply(batches, function(x) x$per_cell_accessible_loci %>% length()))
  accessible_loci_counts_smoothened_list <- split(accessible_loci_counts_smoothened, as.numeric(rep(names(samplesize), samplesize)))
  censored_cor_normalized_variable_mat_list <- lapply(batches, function(x) x$censored_cor_normalized_variable_mat)
  rank_function <- function(i){
    censored_cor_normalized_variable_mat <- censored_cor_normalized_variable_mat_list[[i]]
    accessible_loci_counts_smoothened <- accessible_loci_counts_smoothened_list[[i]]
    nnls_res <- nnls::nnls(censored_cor_normalized_variable_mat,accessible_loci_counts_smoothened)
    nx1_mat <- censored_cor_normalized_variable_mat %*% nnls_res$x
    stepcount = 1
    delta = 1
    extension_parameter=0.99
    original_nx1_mat <- nx1_mat
    current_nx1_mat <- nx1_mat
    next_nx1_mat <- nx1_mat
    while(stepcount < 50000 & delta >= 1e-9){
      stepcount = stepcount + 1
      next_nx1_mat <- extension_parameter*(censored_cor_normalized_variable_mat%*% current_nx1_mat) + (1-extension_parameter) * original_nx1_mat
      delta <- mean(abs(next_nx1_mat-current_nx1_mat))
      current_nx1_mat <- next_nx1_mat
    }
    epitrace_ranked <- rank(current_nx1_mat)
    epitrace_ranked
  }
  if(parallel==T){
    epitrace_ranked <- parallel::mclapply(1:length(censored_cor_normalized_variable_mat_list), mc.cores = min(chunk,ncores),function(i) {
      rank_function(i)
    }) %>% unlist()
  }else{
    epitrace_ranked <- lapply(1:length(censored_cor_normalized_variable_mat_list), function(i) {
      rank_function(i)
    }) %>% unlist()
  }
  final_cells <- colnames(mat2)
  epitrace_norm <- (epitrace_ranked-min(epitrace_ranked))/(max(epitrace_ranked)-min(epitrace_ranked))
  data.frame(
    cell=final_cells,
    EpiTrace = epitrace_norm,
    Accessible_Loci = per_cell_accessible_loci,
    Accessible_Loci_Smooth = accessible_loci_counts_smoothened
  ) -> returndf
  return(returndf)
}






#' AssociationOfPeaksToAge: function for computing each peaks' linear correlation to EpiTrace age.
#' @title AssociationOfPeaksToAge
#'
#' @description function for computing each peaks' linear correlation to EpiTrace age.
#'
#' @details AssociationOfPeaksToAge(epitrace_object,ref_ClockDML_name='AllClock',epitrace_age_name='EpiTraceAge_all',subset_of_cells=NULL,epitrace_age_vector=NULL)
#' @details This is not going to work if the ref_ClockDML_name does not corresponds to epitrace_age_name.
#'
#' @param epitrace_object a seurat object with peaks overlapping reference ClockDML, and computed epitrace age.
#' @param peakSetName a peakset name as in 'assays' of the object. Default to 'peaks'.
#' @param epitrace_age_name name of output epitrace age
#' @param subset_of_cells subset of cell types (Idents) that you want to use in the computation
#' @param epitrace_age_vector alternatively you can have a vector of ages that you would like to test
#' @param parallele switch for parallel						      
#'
#'
#' @return a data frame with locus (the peak), correlation_of_EpiTraceAge, scaled_correlation_of_EpiTraceAge.
#' @export
#' @examples

AssociationOfPeaksToAge <- function(epitrace_object,peakSetName='peaks',epitrace_age_name='EpiTraceAge_all',subset_of_cells=NULL,epitrace_age_vector=NULL,parallel=F){
  if(!is.null(subset_of_cells)){
    subset(epitrace_object,celltype %in% subset_of_cells) -> epitrace_object
  }else{
    epitrace_object -> epitrace_object
  }
  if(DefaultAssay(epitrace_object)!=peakSetName){
    DefaultAssay(epitrace_object) <- peakSetName
  }
  Seurat::GetAssayData(epitrace_object,slot='data') -> peaks_PT_dat
  if(is.null(epitrace_age_vector)){
    epitrace_age_vector <- epitrace_object@meta.data[,epitrace_age_name] %>% as.numeric()
  }
  cell_res_single <- epitrace_object@meta.data %>% as.data.frame()
  # add parallel and avoid explosion
  if(is.null(names(epitrace_age_vector))){
	names(epitrace_age_vector) <- rownames(epitrace_object@meta.data)
  }
  epitrace_age_vector -> age_vec
  is.na(age_vec) -> remove_vec
  peaks_PT_dat -> to_be_associated_mtx
  if(sum(remove_vec)>0){
	age_vec <- age_vec[!remove_vec]
	to_be_associated_mtx <- to_be_associated_mtx[,!remove_vec]
  }  
  if(parallel==F){
	  lapply(c(1:ceiling(dim(to_be_associated_mtx)[1]/1000)),function(x){
	    WGCNA::cor(x = t(to_be_associated_mtx[((1000*(x-1))+1):min(dim(to_be_associated_mtx)[1],1000*x),]), y = age_vec) 
	  }) -> res_list_cor
  }else{
	  parallel::mclapply(c(1:ceiling(dim(to_be_associated_mtx)[1]/1000)),mc.cores = 12,function(x){
	    WGCNA::cor(x = t(to_be_associated_mtx[((1000*(x-1))+1):min(dim(to_be_associated_mtx)[1],1000*x),]), y = age_vec) 
	  }) -> res_list_cor
  }
  res_list_cor %>% unlist -> cor_res_PT

  # cor_res_PT <- WGCNA::cor(x=t(peaks_PT_dat),y=epitrace_age_vector)
  scale(cor_res_PT) -> scaled_cor_res_PT
  names(cor_res_PT) <- rownames(peaks_PT_dat)
  data.frame('locus'=names(cor_res_PT),correlation_of_EpiTraceAge=cor_res_PT,scaled_correlation_of_EpiTraceAge=scaled_cor_res_PT) -> returndf
  return(returndf)
}




#' EpiTraceAge_Convergence: wrapper function for iterative optimization for EpiTrace age (on a specific system)
#' @title EpiTraceAge_Convergence
#'
#' @description wrapper function for iterative optimization for EpiTrace Age. 
#'     The function goes through building an EpiTrace object with given clock reference (the DML). 
#'     Initial derivation of sample age was performed with clock reference only.
#'     Age-association of ChrAcc peaks were performed to extract informative peaks. 
#'     Iterative inference of cell age was performed with an updated clock reference containing 
#'     the top peaks associating with age. 
#'     The interation goes up to given times or when age derivation converge under given error limit. 
#'     Please use **all** peaks instead of **partial** peaks for this optimization. 
#'     Generally this is a maximal extending approach to include all possible peaks correlated with age. 						      
#'
#' @details EpiTraceAge_Convergence (peakSet,matrix,celltype=NULL,min.cutoff=50,lsi_dim=2:50,fn.k.param=21,ref_genome='hg38',sep_string= c(":", "-"),clock_gr_list=clock_gr_list,non_standard_clock=F,qualnum = 10,Z_cutoff=3,mean_error_limit=1e-2)
#' @details Note: SC/Bulk could not be run at the same time as the normalization/filtering process kills single cells in the face of bulk data
#'
#' @param matrix input count matrix of ATAC data, rows=GRanges and cols=samples/single cells.
#' @param peakSet a GenomicRanges object corresponding to the rows of count matrix.
#' @param min.cutoff minimal cutoff for Signac variableFeature calling
#' @param lsi_dim the dimensionalities used in LSI for Signac
#' @param fn.k.param the k parameter used in FindNeighbours
#' @param ref_genome hg38 or hg19. Currently only support these two.
#' @param sep_string the separation string for row names in input matrix to generate ranges. for example, 'chr1:1-2' is c(':','-')
#' @param clock_gr_list the clockDML set, is a list of reference GRanges, each corresponds to a set of clock-like DML or DMR
#' @param non_standard_clock whether is not using non-standard reference clockDML sets.
#' @param qualnum minimal peak/cell number for classifying a 'qualified' cell/peak, default set to 10, can be 1 for including most cells/peaks. 
#' @param Z_cutoff a cutoff for scaled (Z) correlation between peak and cell age. usually set to >= 2.5
#' @param mean_error_limit a limit for the difference between cell ages derived from previous and current iteration. usually set to 1e-2 - 1e-6. 
#' @param ncore_lim limit for using parallel cores in age-peak association, default 12 (set to 1 if no parallel is wanted)
#' @param parallel whether use parallel multicore. In Rstudio, switch off to avoid problems
#' @param iterative_time number of iterations. Usually >10 should be fine. 					      
#'
#' @return a seurat object with EpiTraceAge, Accessibility, and AccessibilitySmoothed. '_iterative' are iteration scores, and '_Clock_initial' are initial mitosis clock scores. 
#' @export
#' @examples



EpiTraceAge_Convergence <- function (peakSet, matrix, celltype = NULL, min.cutoff = 50, lsi_dim = 2:50, fn.k.param = 21, ref_genome = "hg38", sep_string = c(":", "-"), clock_gr = plyranges::reduce_ranges(c(clock_gr_list[[1]],clock_gr_list[[2]])), non_standard_clock = F, qualnum = 10,Z_cutoff=3,mean_error_limit=1e-2,ncore_lim=12,parallel=T,iterative_time=2) {
  
  # debug
  
  # peakSet <- initiated_peaks
  # matrix <- initiated_mm
  # celltype <- NULL
  # ref_genome = 'hg38'
  # non_standard_clock = T
  # mean_error_limit=1e-2
  # min.cutoff = 50
  # lsi_dim = 2:50
  # fn.k.param = 21
  # sep_string = c(":", "-")
  # clock_gr = plyranges::reduce_ranges(c(clock_gr_list[[1]],clock_gr_list[[2]]))
  # non_standard_clock = F
  # qualnum = 1
  # Z_cutoff=3
  # mean_error_limit=1e-2
  # ncore_lim=46
  
  
  # initiate
  original_clk_peakset <- clock_gr
  if (ref_genome  %in% "hg38") {
    original_clk_peakset <- easyLift::easyLiftOver(original_clk_peakset,'hg19_hg38')
  }
  if (ref_genome != "hg19" & ref_genome != "hg38") {
    message("please make double sure your ref genome, peak set and cells are similar.")
  }
  iterative_count = 1
  na_vector_current <- ncol(matrix)
  na_vector_previous <- ncol(matrix)
  error <- rep(1,ncol(matrix))
  mean_error <- mean(error)
  iterative_GR_list <- list('iterative'=original_clk_peakset)
  
  findOverlaps(peakSet,original_clk_peakset)@from %>% unique() -> overlap_with_clk
  matrix[overlap_with_clk,] -> initial_matrix_clk
  peakSet[overlap_with_clk,] -> initial_peakSet_clk
  message('Preparing obj...')
  epitrace_obj <- EpiTrace_prepare_object(initial_peakSet_clk,initial_matrix_clk,celltype,ref_genome = 'hg19',non_standard_clock = T,clock_gr_list = iterative_GR_list,sep_string=sep_string,fn.k.param = fn.k.param,lsi_dim = lsi_dim,qualnum = qualnum,min.cutoff=min.cutoff,run_reduction=F)  # note here I do not switch the ref genome in prepare_obj, for simplicity. 
  message('Finished prepare obj')
  epitrace_obj_original_metadata <- epitrace_obj@meta.data
  message('Estimating age for initialization...')
  if(parallel==F){
    epitrace_obj_age_estimated <- RunEpiTraceAge(epitrace_obj) %>% suppressMessages() %>% suppressWarnings()
  }else{
    epitrace_obj_age_estimated <- RunEpiTraceAge(epitrace_obj,ncores=ncore_lim,parallel=T) %>% suppressMessages() %>% suppressWarnings()
  }
  message('Finished age estimation')
  age_current <- epitrace_obj_age_estimated@meta.data$EpiTraceAge_iterative
  na_vector_current <- is.na(age_current)
  # matrix %>% as.matrix() -> to_be_associated_mtx
  matrix -> to_be_associated_mtx
  next_peakset <- iterative_GR_list$iterative
  
  # res_list <- list() # debug
  while((sum(na_vector_current)<=sum(na_vector_previous)) & (mean_error >= mean_error_limit) & (iterative_count <= iterative_time)){
    message('Iterating ',iterative_count)
    age_previous <- age_current
    na_vector_previous <- na_vector_current
    iterative_count = iterative_count + 1
    remove_vec <- is.na(age_current)
    if(sum(remove_vec)>0){
      age_current <- age_current[!remove_vec]
      to_be_associated_mtx <- to_be_associated_mtx[,!remove_vec]
    }  
    if(parallel==T){
      parallel::mclapply(c(1:ceiling(dim(to_be_associated_mtx)[1]/1000)),mc.cores = ncore_lim,function(x){
        target_mtx <- t(to_be_associated_mtx[((1000*(x-1))+1):min(dim(to_be_associated_mtx)[1],1000*x),]) %>% as.matrix() # otherwise got large cholmod
        WGCNA::cor(x = target_mtx, y = age_current) 
      }) -> res_list_cor
      res_list_cor %>% unlist -> cor_res_PT
    }else{
      lapply(c(1:ceiling(dim(to_be_associated_mtx)[1]/1000)),function(x){
        target_mtx <- t(to_be_associated_mtx[((1000*(x-1))+1):min(dim(to_be_associated_mtx)[1],1000*x),]) %>% as.matrix() # otherwise got large cholmod
        WGCNA::cor(x = target_mtx, y = age_current) 
      }) -> res_list_cor
      res_list_cor %>% unlist -> cor_res_PT
    }
    cor_res_PT[is.na(cor_res_PT)] <- 0 # remove bug 
    scale(cor_res_PT) -> scaled_cor_res_PT
    updated_peakset <- peakSet
    updated_peakset$correlation_of_EpiTraceAge <- cor_res_PT
    updated_peakset$scaled_correlation_of_EpiTraceAge <- scaled_cor_res_PT
    findOverlaps(updated_peakset,original_clk_peakset)@from %>% unique -> peaks_overlap_with_clock
    updated_peakset$locus_type <- 'Others'
    updated_peakset$locus_type[peaks_overlap_with_clock] <- 'Chronology'
    list('iterative'=updated_peakset[(updated_peakset$scaled_correlation_of_EpiTraceAge < -1 * Z_cutoff | updated_peakset$locus_type %ni% 'Others') %>% as.vector() | updated_peakset$peakId %in% next_peakset$peakId,]) -> iterative_clock_gr_list # update new age-associated peaks correlated to mitosis score
    findOverlaps(peakSet,iterative_clock_gr_list[[1]])@from %>% unique() -> overlap_with_clk
    matrix[overlap_with_clk,] -> initial_matrix_clk
    peakSet[overlap_with_clk,] -> initial_peakSet_clk
    # epitrace_obj_iterative <- EpiTrace_prepare_object(initial_peakSet_clk,initial_matrix_clk,celltype,ref_genome = 'hg19',non_standard_clock = T,clock_gr_list = iterative_clock_gr_list) %>%suppressMessages() %>% suppressWarnings()
    if(parallel==T){
      if(nrow(initial_matrix_clk)*ncol(initial_matrix_clk)>50000*2000){
        batch_size=1000
      }else{
        batch_size=400 # size=2000,batch=NULL,subsamplesize=2000
      }
      EpiTraceAge(initial_matrix_clk,parallel=parallel,ncores=ncore_lim,size=batch_size,subsamplesize=batch_size) %>%suppressMessages() %>% suppressWarnings() -> res1
    }else{
      EpiTraceAge(initial_matrix_clk,parallel=F,ncores=1) %>%suppressMessages() %>% suppressWarnings() -> res1
    }
    
    colnames(res1)[2:4] <- c('EpiTraceAge_iterative','Accessibility_iterative','AccessibilitySmooth_iterative')
    epitrace_obj_original_metadata_update <- left_join(epitrace_obj_original_metadata,res1)
    epitrace_obj_iterative_age_estimated <- epitrace_obj_age_estimated
    rownames(epitrace_obj_original_metadata_update) <- rownames(epitrace_obj_original_metadata)
    epitrace_obj_iterative_age_estimated@meta.data <- epitrace_obj_original_metadata_update
    # 
    # epitrace_obj_iterative_age_estimated <- RunEpiTraceAge(epitrace_obj_iterative)%>%suppressMessages() %>% suppressWarnings()
    age_current <- epitrace_obj_iterative_age_estimated$EpiTraceAge_iterative
    error <- (age_current - age_previous)
    tryCatch({
      error[is.infinite(error)] <- 0
    },error=function(e){message('')})
    tryCatch({
      error[is.na(error)] <- 1
    },error=function(e){message('')})
    mean_error = mean(abs(error),na.rm=T)
    na_vector_current <- is.na(age_current)
    age_current -> age_previous 
    epitrace_obj_iterative_age_estimated@misc$iterative_count <- iterative_count
    epitrace_obj_iterative_age_estimated@misc$mean_error <- mean_error
    message('mean_error = ',mean_error)
    next_peakset <- iterative_clock_gr_list$iterative
    # res_list[[iterative_count]] <- list('iterative_clock_gr_list'=iterative_clock_gr_list,mean_error=mean_error,age_current=age_current,iterative_count=iterative_count) # debug
  }
  epitrace_obj_iterative_age_estimated$EpiTraceAge_Clock_initial <- epitrace_obj_age_estimated@meta.data$EpiTraceAge_iterative
  epitrace_obj_iterative_age_estimated$Accessibility_initial <- epitrace_obj_age_estimated@meta.data$Accessibility_iterative
  epitrace_obj_iterative_age_estimated$AccessibilitySmooth_initial <- epitrace_obj_age_estimated@meta.data$AccessibilitySmooth_iterative
  # return(list(epitrace_obj_iterative_age_estimated=epitrace_obj_iterative_age_estimated,res_list=res_list)) # debug
  return(epitrace_obj_iterative_age_estimated)
}
  
  









