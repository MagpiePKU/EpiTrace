

RunEpiTracePhylogeny
--------------------

Input: 

- epitrace_object : an EpiTrace object 
- min.cutoff: used in Signac pre-processing step. Default to 50. 
- run_reduction: whether run Signac dimensionality reduction. Default to TRUE. 

Output: a list consisting (1) the BuildTree object; (2) the phylogenetic tree from ape; and (3) a ggtree object. 

The function is a wrapper of Seurat build-tree utility. 

**We plan to depreciate the ggtree utility soon** 

  