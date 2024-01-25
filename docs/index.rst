

EpiTrace - Estimating cell age with bulk and single-cell ATAC-seq
============================================================

.. image:: https://github.com/MagpiePKU/EpiTrace/assets/7695551/b1104b09-bbbb-420f-bac0-8241e828e2bc
   :width: 300px
   :align: left



**EpiTrace** is an R package for estimating cell age from single-cell ATAC-seq data. It takes an approximation approach to infer the relative mitosis (replicative) age -- the number of mitosis a cell has undergone. It does so by measuring the total opened reference "clock-like" genomic loci. On these loci, heterogeneity of chromatin accessibility decreases as the cell ages. The chromatin accessibility-based mitosis age inferred by EpiTrace adds a time domain to the single-cell sequencing data. It complements somatic mutation, RNA velocity and stemness predictions to predict the cell evolution trajectory with improved precision and power :cite:p:`Xiao2022`. 




EpiTrace's key applications
^^^^^^^^^^^^^^^^^^^^^^^^^^^
- estimate the mitosis age of single cell or bulk sample.  
- identify age-dependent biological events including shifts in chromatin accessibility, transcription factor activity, and transcriptomic changes. 
- timing developmental and mutational events. 
- infer the total ATAC or ChIP-seq signal over a given reference genomic region set. 


Citing EpiTrace
^^^^^^^^^^^^^^^

If you include or rely on EpiTrace when publishing research, please adhere to the
following citation guide:

**EpiTrace algorithm and the ClockDML**

If you use the *algorithm* and/or *ClockDML* (*clock-like differential methylated loci*), including cross-species lift-over clock-like loci generated from the ClockDML, cite

.. code-block:: bibtex

    @article {Xiao2022.05.12.491736,
	author = {Yu Xiao and Wan Jin and Lingao Ju and Jie Fu and Gang Wang and Mengxue Yu and Fangjin Chen and Kaiyu Qian and Xinghuan Wang and Yi Zhang},
	title = {Tracking single cell evolution via clock-like chromatin accessibility},
	elocation-id = {2022.05.12.491736},
	year = {2024},
	doi = {10.1101/2022.05.12.491736},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/01/06/2022.05.12.491736},
	eprint = {https://www.biorxiv.org/content/early/2024/01/06/2022.05.12.491736.full.pdf},
	journal = {bioRxiv}}


**G-quadruplex region**

If you use the *G-quadruplex region* for cell/sample age estimation, cite

.. code-block:: bibtex

    @article {Jin2024.01.06.574476,
	author = {Wan Jin and Jing Zheng and Yu Xiao and Lingao Ju and Fangjin Chen and Jie Fu and Hui Jiang and Yi Zhang},
	title = {A universal molecular mechanism driving aging},
	elocation-id = {2024.01.06.574476},
	year = {2024},
	doi = {10.1101/2024.01.06.574476},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/01/06/2024.01.06.574476},
	eprint = {https://www.biorxiv.org/content/early/2024/01/06/2024.01.06.574476.full.pdf},
	journal = {bioRxiv}}



Support
^^^^^^^
Have a question or would like to start a new discussion? Found a bug or would like to see a feature implemented? Feel free to submit an
`issue <https://github.com/MagpiePKU/EpiTrace/issues/new>`_.
Your help to improve EpiTrace is highly appreciated. 

Planned updates
^^^^^^^^^^^^^^^
Tutorial on (1) using G-quadruplex as reference clock-like loci, and (2) timing mutational event during oncogenesis.  


.. toctree::
   :caption: Main
   :maxdepth: 1
   :hidden:

   about
   installation
   release_notes
   references


.. toctree::
   :caption: Tutorials
   :maxdepth: 1
   :hidden:

   Bulk_ATAC
   scATAC_cIPSC
   scATAC_drosophila
   mtscATAC_CD34


.. toctree::
   :caption: Functions
   :maxdepth: 1
   :hidden:

   Init_Peakset 
   Init_Matrix
   EpiTrace_prepare_object
   RunEpiTraceAge
   EpiTraceAge_Convergence
   AssociationOfPeaksToAge
   RunEpiTracePhylogeny





.. |br| raw:: html

  <br/>

.. |dim| raw:: html

   <span class="__dimensions_badge_embed__" data-id="pub.1129830274" data-style="small_rectangle"></span>
   <script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>
