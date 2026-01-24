Interpreting EpiTrace Age
================================================================================

Understanding Bulk vs Single-Cell EpiTrace Age
================================================================================

.. important::
   **EpiTrace age interpretation differs between bulk and single-cell ATAC-seq data**

Key Points
--------------------------------------------------------------------------------

* **Single-cell ATAC-seq**: EpiTrace age is **positively correlated** with mitotic age
  - Higher EpiTrace age = Older cells (more mitotic divisions)
  - No transformation needed

* **Bulk ATAC-seq**: EpiTrace age is **negatively correlated** with mitotic age
  - Higher EpiTrace age = Younger cell populations
  - Transformation: ``1 - rank/max(rank)`` may be applied for visualization

Why the Difference?
================================================================================

The difference arises from how chromatin accessibility patterns are aggregated:

1. **Single-cell data**: EpiTrace measures accessibility per cell
   - Cells with more mitotic divisions have lost clock-like accessibility
   - The algorithm directly measures this loss

2. **Bulk data**: EpiTrace measures average accessibility across cell populations
   - Populations with more mitotic cells have lost clock-like accessibility
   - The inverse relationship emerges from population averaging

How EpiTrace Functions Handle This
================================================================================

EpiTrace functions do **NOT** automatically detect or transform based on data type:

* ``RunEpiTraceAge()`` - Computes age without considering bulk vs single-cell
* ``EpiTraceAge_Convergence()`` - Iterative refinement, no automatic transformation
* ``EpiTraceAge()`` - Core age calculation function

**You must manually apply the appropriate interpretation for your data type.**

Practical Examples
================================================================================

Single-Cell Data Example
--------------------------------------------------------------------------------

.. code-block:: r

    # Load single-cell ATAC data
    data("scATAC_example")

    # Run EpiTrace
    obj <- EpiTrace_prepare_object(peakset, matrix, celltype)
    obj <- RunEpiTraceAge(obj)

    # Interpret: Higher EpiTraceAge = Older cells
    # No transformation needed
    plot(obj$EpiTraceAge_AllClock, main="EpiTrace Age (higher = older)")

Bulk Data Example
--------------------------------------------------------------------------------

.. code-block:: r

    # Load bulk ATAC data
    data("bulkATAC_example")

    # Run EpiTrace
    obj <- EpiTrace_prepare_object(peakset, matrix, celltype)
    obj <- RunEpiTraceAge(obj)

    # Interpret: Higher EpiTraceAge = Younger population
    # For visualization, you may want to invert:
    obj$EpiTraceAge_Inverted <- 1 - obj$EpiTraceAge_AllClock
    plot(obj$EpiTraceAge_Inverted, main="EpiTrace Age (higher = older)")

When to Use Transformation
================================================================================

**Do NOT transform** (use raw EpiTrace age):
* Single-cell data analysis
* Comparing relative age ordering within cells
* Using age as a continuous variable in downstream analysis

**Consider transforming** (e.g., ``1 - rank/max(rank)``):
* Bulk data visualization
* Making bulk results intuitive ("higher = older")
* Publication figures where biological age increases from left to right

Common Pitfalls
================================================================================

1. **Mixing bulk and single-cell data**: Always analyze them separately, then compare results

2. **Forgetting to document**: Always note whether you're using raw or transformed EpiTrace age

3. **Comparing across platforms**: Bulk and single-cell EpiTrace ages are not directly comparable

References
================================================================================

For more details, see:
* Xiao et al. (2024) Nature Biotechnology - Methods section
* `Bulk_ATAC` tutorial - Bulk data examples
* `scATAC_cIPSC` tutorial - Single-cell data examples
