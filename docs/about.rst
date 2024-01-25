About EpiTrace
--------------

Single cell chromatin accessibility sequencing (scATAC) reconstructs cell developmental trajectory by phenotypic similarity. However, inferring the exact developmental trajectory is challenging. While excellent tools such as RNA velocity, stemness prediction and metabolic labeling exist for determining the cell evolution trajectory on the manifold of phenotypes – usually described as Waddington’s landscape – for single cell RNA-seq datasets, no comparable methods exist for scATAC. State-of-the-art, similarity-based lineage deduction methods would be limited when phenotypes are fluidic, such as in dedifferentiation or oncogenesis. On the other hand, mutation-based lineage tracing methods, for example, using mitochondrial SNPs, which track the phylogeny of cells over divisions, are highly accurate, yet their temporal resolution is restrained by the low natural mutation rate.

Introducing a time domain into these scATAC data would enable studying dynamic processes during development and aging. Such information, if exists, would likely to be intimately related cell replication :cite:p:`Hayflick1965`. Historically, epigenetic modification such as DNA methylation has been shown to be associated with biological age :cite:p:`Horvath2013` as well as mitosis :cite:p:`Youn2018`. In our work, we have identified that chromatin accessibility on certain genomic regions also exhibit clock-like behavior, such that the heterogeneity of chromatin accessibility on these regions reduces across mitosis. EpiTrace, described in :cite:p:`Xiao2022`, leverages this phenomenon to estimate single cell ages from a given set of reference clock-like region. 


The replication-associated chromatin accessibility change
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: https://github.com/MagpiePKU/EpiTrace/assets/7695551/8bdca538-7530-4535-8193-e6d517ae6c81
   :width: 250px
   :align: left

Theoretically, the trajectory of epigenome reprogramming during cell state transformation, either in aging or development, could be decomposed into two orthogonal processes: cell differentiation (or specification) and cell replication. While the cell differentiation components might be orthogonal (*distinct*) between cell types or organisms, there must be a similar cell replication component. Otherwise, there would not be a common aging-associated epigenomic biomarker, which is contradictory to our and other groups’ findings. 

We have determined that DNA G-quadruplex stimulates both RNA transcription and DNA replication, inducing local transcription-replication interaction which delays genome replication. The delayed genome replication impairs epigenomic modification transmission from parental to the newly synthesized daughter genome, resulting in *DNA hypomethylation* and loss of heterochromatin. This creates a more permissive environment for G-quadruplex formation in the subsequent generation. As a result, DNA G-quadruplexes gradually accumulate across mitosis to erode the local epigenome. Accompanying with this, chromatin around putative G-quadruplex sequences gradually *opens* across replication. Please refer to :cite:p:`Jin2024` for the details on this molecular mechanism.  


The phenotypic-neutral, reference clock-like loci
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The reference clock-like loci, as described earlier, is a set of genomic regions that shows clock-like chromatin accessibility. Namely, they undergo *irreversible* chromatin accessibility changes across mitosis. In other words, on a population scale, they either monotonically open or close during aging. Therefore, regardless of their initial state, the overall chromatin accessibility profile on them eventually converge into a defined state. Examples of such loci includes (1) the clock-like differentially methylated loci in human and mouse; (2) homologous region of clock-like differentially methylated loci in other species; and (3) putative G-quadruplex sequences. 

On the single cell scale, reduction of global heterogeneity of chromatin accessibility across a set of reference genomic loci could be approximated by the fraction of opened reference loci :cite:p:`Xiao2022`. EpiTrace estimates this index by using HMM-mediated smoothing on a cell-to-cell similarity matrix (conceptually similar to :cite:p:`Gulati2020`). 


Identifying cell-type-specific clock-like loci
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a particular cell lineage, the epigenomic shift during aging could be cell-type dependent. Furthermore, between different cell clones, the epigenomic shift on similar putative clock-like loci could be stochastic. Nevertheless, we can still infer the cell-type-specific clock-like loci by correlating their chromatin accessibility to the approximated cell age. EpiTrace does so by: 

- compute the correlation of chromatin accessibility of each loci to single cell age estimation (initially estimated by using the reference clock-like loci); 
- search for loci that shows best correlation; 
- combine these loci with the previous reference clock-like loci to build a new set of reference; 
- estimate single cell age using the new set of reference loci. 

.. image:: https://github.com/MagpiePKU/EpiTrace/assets/7695551/bf9c031a-e9ee-4714-bf47-b3a309935c66
   :width: 600px
   :align: center

This process is iterated to update the reference loci until the estimated single cell age converges. In practice, we find ~10 iterations usually is sufficient to give satisfactory result. 

See :cite:p:`Xiao2022` for a detailed exposition of the methods.
