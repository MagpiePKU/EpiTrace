# EpiTrace
 Computing cell age with bulk and single-cell ATAC-seq data   

 EpiTrace takes an approximation approach to infer single cell age from single cell ATAC data. It infers single cell age by measuring the total opened reference genomic loci. On these loci, heterogeneity of chromatin accessibility decreases as the cell ages.   

 EpiTrace firstly algorithmically determine a set of tool genomic loci, on which the total chromatin accessibility (reads) shows maximal correlation to the total opened reference genomic loci. Then, the total chromatin accessibility on this set of tool genomic loci are used as an intermediate tool variable to approximate cell age.  
 
 EpiTrace documentation is now on `readthedocs`. 
 
 For descriptions, function references, and tutorials, visit https://epitrace.readthedocs.io 

 Maintainer: Zhang Yi <c.sinensis@gmail.com>      

### Installation
```
if(!require(pak)){
    install.packages("pak")
}
library(pak)
pak::pkg_install('caleblareau/easyLift')
pak::pkg_install('MagpiePKU/EpiTrace')  
```

