Clock Reference Options
================================================================================

EpiTrace provides multiple reference clock-like loci sets for age estimation.

Available Clock Sets
================================================================================

Standard Clock Sets (included with the package)
--------------------------------------------------------------------------------

1. **Mitosis** (616 regions)
   - Contains only **hypermethylation** sites
   - Shows strong positive correlation with mitotic age
   - Best for: Detecting recent cell divisions, identifying highly proliferative cells
   - Use when: Focus on mitosis-specific age estimation

2. **Chronology** (125,625 regions)
   - Contains clock-like DML/DMR from multiple studies
   - Excludes regions overlapping with Mitosis clock
   - Mix of hyper- and hypomethylation patterns
   - Best for: General age estimation across cell types
   - Use when: Balanced view of cellular aging

3. **solo_WCGW** (1,669,720 regions)
   - Solitary WCGW (WCGW = Ws/C/Gs) motifs in partially methylated domains
   - Contains only **hypomethylation** sites
   - Large number of regions (1.6M vs 616 for Mitosis)
   - Best for: Research on WCGW-specific patterns, comparative studies
   - Use when: Investigating hypomethylation-specific age signatures

4. **AllClock** (Mitosis + Chronology)
   - Combined set of Mitosis and Chronology clocks (126,241 regions)
   - Default when using standard clock with ``RunEpiTraceAge()``
   - Best for: Most applications, comprehensive age estimation
   - Use when: Standard EpiTrace analysis

Which Clock Should You Use?
================================================================================

For Most Applications (Recommended)
--------------------------------------------------------------------------------

Use **AllClock** (the default):

.. code-block:: r

    # Default: uses Mitosis + Chronology
    obj <- EpiTrace_prepare_object(peakset, matrix, celltype,
                                   clock_gr_list = clock_gr_list,
                                   standard_clock = TRUE)
    obj <- RunEpiTraceAge(obj)

For Mitosis-Specific Age Estimation
--------------------------------------------------------------------------------

Use **Mitosis** clock only:

.. code-block:: r

    # Use only Mitosis clock
    custom_clock <- list(MitosisClock = clock_gr_list[['Mitosis']])
    obj <- EpiTrace_prepare_object(peakset, matrix, celltype,
                                   clock_gr_list = custom_clock,
                                   non_standard_clock = TRUE)

For Hypomethylation Studies
--------------------------------------------------------------------------------

Use **solo_WCGW** clock:

.. code-block:: r

    # Use solo_WCGW for hypomethylation-focused analysis
    custom_clock <- list(solo_WCGW = clock_gr_list[['solo_WCGW']])
    obj <- EpiTrace_prepare_object(peakset, matrix, celltype,
                                   clock_gr_list = custom_clock,
                                   non_standard_clock = TRUE)

Why is solo_WCGW Not Used in Tutorials?
================================================================================

The tutorials use **AllClock** (Mitosis + Chronology) rather than solo_WCGW for several reasons:

1. **Computational efficiency**
   - solo_WCGW has 1.6M regions vs 616 for Mitosis
   - Processing time increases significantly with more regions
   - Most applications don't need this granularity

2. **Biological interpretation**
   - Mitosis + Chronology provides a balanced view
   - Combines hyper- and hypomethylation patterns
   - Validated in multiple cell types and systems

3. **Validation**
   - AllClock (Mitosis + Chronology) is extensively validated
   - Shown to work well across different datasets
   - solo_WCGW is primarily for specialized research questions

When to Consider solo_WCGW
================================================================================

Use solo_WCGW if:

* Studying WCGW motif-specific aging patterns
* Investigating hypomethylation-driven age changes
* Comparing hyper- vs hypomethylation age signatures
* Working with model systems where WCGW patterns are prominent
* Exploring age-associated chromatin changes in partially methylated domains

Example Comparative Analysis
================================================================================

.. code-block:: r

    # Compare different clock sets
    clocks_to_compare <- list(
      Mitosis = clock_gr_list[['Mitosis']],
      Chronology = clock_gr_list[['Chronology']],
      solo_WCGW = clock_gr_list[['solo_WCGW']]
    )

    results <- list()
    for (clock_name in names(clocks_to_compare)) {
      obj <- EpiTrace_prepare_object(peakset, matrix, celltype,
                                     clock_gr_list = list(clock = clocks_to_compare[[clock_name]]),
                                     non_standard_clock = TRUE)
      obj <- RunEpiTraceAge(obj)
      results[[clock_name]] <- obj$EpiTraceAge_clock
    }

    # Compare correlations
    cor(results[['Mitosis']], results[['Chronology']])
    cor(results[['Mitosis']], results[['solo_WCGW']])

References
================================================================================

For more details on clock set construction:
* Xiao et al. (2024) Nature Biotechnology - Methods section
* Yang et al. (2016) - Mitosis clock DML
* Youn and Wang (2018) - Chronology clock DML
* Zhou et al. (2018) - solo_WCGW in partially methylated domains
