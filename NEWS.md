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
