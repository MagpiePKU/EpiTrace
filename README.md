# EpiTrace
 Computing cell age with bulk and single-cell ATAC-seq data   

 EpiTrace takes an approximation approach to infer single cell age from single cell ATAC data. It infers single cell age by measuring the total opened reference genomic loci. On these loci, heterogeneity of chromatin accessibility decreases as the cell ages.   

 EpiTrace firstly algorithmically determine a set of tool genomic loci, on which the total chromatin accessibility (reads) shows maximal correlation to the total opened reference genomic loci. Then, the total chromatin accessibility on this set of tool genomic loci are used as an intermediate tool variable to approximate cell age.  
 
 For example code, see vignettes/Hematopoiesis_2016_demo.notebook.Rmd or vignettes/Hematopoiesis_2016_demo.notebook.nb.html  

 Maintainer: Zhang Yi <c.sinensis@gmail.com>      

### Installation
```
library(devtools)   
devtools::install_github('MagpiePKU/EpiTrace')    
```

### User's handbook
EpiTrace tutorial is now online! Please visit https://epitrace.readthedocs.io 
