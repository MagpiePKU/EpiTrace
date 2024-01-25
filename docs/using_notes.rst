Known issues and workarounds
----------------------------

Skewed dataset 
''''''''''''''
Please note that EpiTrace takes an approximation approach to infer single cell age. In its current form, the estimation of tool genomic loci during first step of EpiTrace algorithm might be affected by the data. Particularly, when single cell composition is biased in the dataset, for example, to have 99% of one type of cell and only 1% of others, then the measurement could be incorrect. We encourage users to use balanced dataset for cell age estimation, for example down-sampling the over-represented cells or simply avoid using the overrepresented sample.  
  

Low sequencing depth 
''''''''''''''''''''
Low sequencing depth results in higher variation in age estimation. Currently, there is no good workaround over this problem. However, this problem should not be specifically affecting EpiTrace performance as these data are generally of less usage in other applications as well. 


