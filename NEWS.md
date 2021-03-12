# PAMES v2.7.1
- fix `select_informative_regions_ext` with flag `return_info=TRUE`
# PAMES v2.7.0
- renamed parameter `platform` to `ref_table`
# PAMES v2.5.0
- added new function (select_informative_sites_ext.R) that:
  - considers also control data distribution
  - use custom platform data (useful to filter out chromosomes X and Y)
  - return named vectors (names are either derived from rownames of tumor_data 
  or created with format "CpG_xxxxxxx")
- compute_purity now requires platform data (not just "450k" or "27k")
  addressing issue[#7](https://github.com/cgplab/PAMES/issues/7))

# PAMES v2.4.0
* deprecated parameter na_threshold in favor of min_sample_frac `compute_AUC` 
* added simplify parameter: if false return a data.frame reporting fraction of NA samples per site

# PAMES v2.3.5
* `compute_AUC` now returns a vector with the names of the CpG probes.

# PAMES v2.3.4
* Add `platform` parameter to `compute_purity` to address issue [#7](https://github.com/cgplab/PAMES/issues/7))

# PAMES v2.3.3
* Fix bug in `reduce_to_regions` (see issue [#6](https://github.com/cgplab/PAMES/issues/6))

# PAMES v2.3.1
* Fix bug in `compute_AUC`

# PAMES v2.3.0
* Add `percentiles` parameter for a more flexible selection of sites

# PAMES v2.2.0
* Even faster cluster reduction
* Fix bug: "top" method retrieved same sites of "hyper"

# PAMES v2.1.0
* Faster and simpler cluster reduction in `select_informative_sites`
* Removed internal functions:
    - `is_too_close`
    - `cluster_reduction`
* Added variable `method` to `select_informative_sites` and `select_informative_regions`
to select **even**, **top**, **hyper** or **hypo** methylated sites.

# PAMES v2.0.0
* PAMES have been rewritten to make the code simpler and the analysis run faster. 
Thanks to [GMFranceschini](https://github.com/GMFranceschini) for suggesting to
use `GenomicRanges`.

* Removed function:
    - `compute_islands_indexes.R`
* Renamed functions:
    - `reduce_to_islands.R` -> `reduce_to_regions.R`
* Renamed data:
    - `cpg_islands_df` -> `cpg_islands`
* Removed data:
    - `bs_toy_indexes`

# PAMES v1.1.0

* BUG correction: setting na_threshold to be greater than 1 didn't generate an error
(and resulted as having na_threshold set to 0.99)

# PAMES v1.0.1

* Faster compute_AUC thanks to Wilcoxon method

# PAMES v1.0.0

## New Features

* Added parallel computation to `compute_AUC` with argument `ncores`.

## Updated dataset

* `cpg_islands` has been renamed as `cpg_islands_df`
* `bs_tumor_toy_data` and `bs_control_toy_data` has been merged into `bs_toy_matrix`

## Changed behaviour

* All beta tables are now required to be percentage values.
