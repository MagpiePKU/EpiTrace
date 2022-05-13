# EpiTrace
 Computing cell age with bulk and single-cell ATAC-seq data   
 
 For example code, see vignettes/Hematopoiesis_2016_demo.notebook.Rmd or vignettes/Hematopoiesis_2016_demo.notebook.nb.html  

 Maintainer: Zhang Yi <synapse@pku.edu.cn>      
 
 


### Updates   
20220510   
1. Iteration code has been updated. Please use EpiTraceAge_Convergence().   


### Citation
Xinghuan Wang, Wan Jin, Gang Wang, Lingao Ju, Fangjin Chen, Kaiyu Qian, Yu Xiao and Yi Zhang, *Tracking single cell evolution via clock-like chromatin accessibility*, (2022) BIORXIV/2022/491736      



![Screenshot](demo.png)
 
### Installation
```
library(devtools)   
devtools::install_github('MagpiePKU/EpiTrace')    
```

### Basic Usage (Example)  
```
library(EpiTrace)  
data("clock_gr_list")  
data("hematopoiesis.2018")   
initiated_peaks <- Init_Peakset(inputpeak)  
inputmatrix <- Init_Matrix(cell_metadata$celltype,initiated_peaks,inputmatrix)
epitrace_obj <- EpiTrace_prepare_object(initiated_peaks,inputmatrix,cell_metadata$celltype,ref_genome = 'hg19',non_standard_clock = F,clock_gr_list = clock_gr_list)  
epitrace_obj_age_estimated <- RunEpiTraceAge(epitrace_obj)  
phylotree_res <- RunEpiTracePhylogeny(epitrace_obj_age_estimated)  
```

### To plot a phylogenetic tree
```
phylotree_res_myeloid <- RunEpiTracePhylogeny(subset(epitrace_obj_age_estimated,celltype %in% c('HSC','MPP','CMP','GMP','Monocyte')))   
mitosis_tree <- phylotree_res_myeloid[['MitosisClock']][[2]]  
mitosis_tree <- ape::root(mitosis_tree,outgroup='HSC')  
plot(mitosis_tree)  
```

### To compute age-peak association  
```
associated_res_HSC <- AssociationOfPeaksToAge(subset(epitrace_obj_age_estimated,celltype %in% 'HSC'))   
```
 

 
 
 
