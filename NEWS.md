# cjbart 0.3.0

* Rename arguments in `het_vimp()` from `model` to `imces` and `outcomes` to `levels` for consistency across functions
* Change `het_vimp()` so it returns a list, including the results and an attribute-level look-up table
* Added ability to run `het_vimp()` calculation in parallel
* Fixed bug when predicting continuous conjoint outcomes
* Add new `pIMCE()` function to calculate population individual-level marginal component effects
* Add new `AMCE()` function that includes Bayesian uncertainty estimate
* Add check for unique levels across attributes. If not, the attribute-level values are modified to ensure uniqueness
* Add functionality to improve quality of outputs
* Update the vignette to reflect new functionality
* Add ability to drop attributes from IMCE analysis
* Remove dependencies not needed

# cjbart 0.2.2

* Re-release after dependencies back on CRAN

# cjbart 0.2.1

* Minor internal changes to keep pace with R development
* Minor changes to vignette wording

# cjbart 0.2.0

* Addresses unbalanced chunk delimiters in package vignette
* Adds new functionality for analyzing heterogeneity after IMCE prediction using random forest variable importance measures

# cjbart 0.1.0

* First release including all core functionality
* Includes wrapper of BART estimation, IMCE prediction, and generics
