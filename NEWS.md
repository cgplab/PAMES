# PAMES v2.2.0
* Even faster cluster reduction
* Fix bug: "top" method retrieved same sites of "hyper"

# PAMES v2.1.0
* Faster and simplier cluster reduction in `select_informative_sites`
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
