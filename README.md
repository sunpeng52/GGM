# GGM
R codes for "Sparse positive-definite estimation for covariance matrices with repeated measurements"

### Basic Functions
`functions.R` & `functions_2.R`: functions to set up covariance matrices, run ADMM, etc.

### Figure 1 & Figure 2

`mm_generate_data.R`: generate simulation data for Model 1 and Model 2 with p =100, 200.

`mm_simulation_100.R` & `mm_simulation_200.R`: apply our proposed procedures to the simulated data.

`mm_rate_100.R`: compute fpr and tpr for Figure 2.

`00 mm_analysis_100.R`: generate Figure 1.

`00 mm_roc_plot.R`: generate Figure 2.

### Figure 3

`mm_comparison_generate_data.R`: generate simulation data for Model 3 and Model 4.

`mm_comparison_simulation.R`: apply our proposed procedures to the simulated data.

`mm_comparison_analysis.R`: compute several metrics for Figure 3.

`00 mm_comparison_plot.R`: generate Figure 3.

### Table 1 in the Supplementary materials

`mm_comparison_soft_admm.R`: simulation for comparison between ADMM and soft thresholding, results are summarized in Table 1 in the supplementary materials.

### Real Data Analysis 

`01_Real_Analysis_New.R`: code for real data analtsis, generating Figure 4.
