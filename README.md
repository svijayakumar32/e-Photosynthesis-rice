### e-Photosynthesis-rice
This repository contains the code and data necessary to reproduce the results presented in the paper by Vijayakumar, S., Wang, Y.,Lin, H.C., Carmo-Silva, E., Long, S.P. and Taylor, S. H. entitled "Tailoring a dynamic model of photosynthetic metabolism towards greater carbon assimilation in rice".

## Overview
This repository provides a generalisable workflow for calibration, parameter optimisation, and input scaling that can be applied to any plant species with available gas exchange and kinetic data.
Our paper describes the steps for optimising resource investment among photosynthetic enzymes to maximise carbon assimilation in a model of C3 photosynthetic metabolism calibrated to rice, summarised in `pipeline_simple_flowchart.png`. 

![pipeline_simple_flowchart](https://github.com/user-attachments/assets/6e3b0a39-c109-4f79-8831-0ecfef8d41b8)


This version of e-Photosynthesis integrates species-specific temperature dependencies for the balance of oxygenation and carboxylation (Vomax/Vcmax) and ribulose-1,6-bisphosphate carboxylase/oxygenase (Rubisco) activation, as well as accounting for tight binding Rubisco inhibitors.
The workflow begins by combining species-specific equations for calculating temperature dependencies of Rubisco catalytic properties with leaf-level gas exchange measurements for Oryza sativa cv. IR64 to derive photosynthetic parameters describing Calvin-Benson-Bassham (CBB) cycle activity, i.e. Vcmax and J. 
These photosynthetic parameters are used to re-scale enzyme activities in e-Photosynthesis before running the model to identify redistributions that are optimal for CO2 assimilation at different [CO2] levels. Following this, additional analyses are applied to evaluate strategies for increasing photosynthesis through overexpressing subsets of 2-6 enzymes to engineer improved photosynthesis under drought stress, current ambient [CO2] and future elevated [CO2].

## Files
The files within the main repository contains scripts, variables and data adapted from the version of e-Photosynthesis stored in https://github.com/cropsinsilico/C3-metabolic-and-leaf-model.
All code licensing information is detailed in `LICENSE`.

Custom scripts and data files are described below:
- `Calculate_Cc.m` - calculates Cc values equivalent to Ca for the FvCB model fit
- `CalculateGrossAssimilation.m` - calculates gross CO2 assimilation
- `Einput7.txt` - a list of the original enzyme activities for e-Photosynthesis model
- `EPS_Drive_GRNs` - calls EPS_Drive to calculate assimilation rates using a specific list of enzyme activities as input.
- `enzyme_adjustment_test_new`- runs sensitivity analysis to determine assimilation rates associated for random fold changes (n = 2000) of various enzyme combinations drawn between FC = 1 (no change) and FC = 1.25 for Rubisco and FC = 2.8 for SBPase, aldolase, PRK, FBPase and TK at Cc = 130 umol mol-1, 250 umol mol-1 or Cc = 360 umol mol-1
- `FvCB_check_Vcmax_adj.m` - checks transition between Rubisco-limited and RuBP limited phases in A/Cc fit
- `Jmax_adj_simple.m` - runs model calibration to find optimal αEnzymes by minimising SSR between assimilation rates of FvCB and e-Photosynthesis models in the RuBP-regeneration limited range of [CO2]
- `MW&Kcat.txt` - a list of molecular weights and kcat values for photosynthetic enzymes
- `PR_constraints.txt` - a list of photorespiratory enzyme activities (optimised Vmax) obtained from running `job_array_gpmain_rice_129.sh`
- `PR_constraints_protein.txt` - a list of photorespiratory enzyme concentrations (optimised protein) converted from `PR_constraints`
- `Vcmax_adj_simple.m` - runs model calibration to find optimal αRubisco by minimising SSR between assimilation rates of FvCB and e-Photosynthesis models in the Rubisco limited range of [CO2]
  
In addition to this, the remaining files are stored in various directories defining their function:

CO_2 response - for fitting leaf-level gas exchange measurements to the Farquhar-vonCaemmerer-Berry (FvCB) model

The Data sub-folder contains:
- `Gas_exchange_measurement_WT_plants.xls` - Excel workbook containing 8 sets of leaf level gas exchange measurements for wild type Oryza sativa cv. IR64
- Filenames beginning with `IR64-A009-07-33-05` contain data for each curve 
 
The Scripts sub-folder contains:
- `msuRACiFit_rice.R` - R script containing fitting procedure for obtaining photosynthetic parameters

Live_Scripts - for running interactive executable versions of scripts - combining code, output, and formatted text in a single MATLAB environment

Outputs - for storing outputs of the enzyme optimisation within subfolders:
- `Data` - contains `Optimization_Procedure.txt`, `Results_optimization_rice_new_3.xlsx`, and `Series_of_job_submissions.txt`
- `Enzymes` - BestMatrix gives the optimal distribution of Vmax values for 67 photosynthetic enzymes (V1-V59)
- `Metabolites` - dplot gives change in metabolite concentrations (which reach steady state at the end of the optimisation)
- `Workspaces` - MATLAB workspaces saved after running gpmain simulation

Plotting - for plotting results of analyses 
- `Farq_ePhoto_comparison_new.m` - plot A/Cc curve for measured average parameters from gas exchange data and using Ac and Aj equations from Vcmax_adj_simple and Jmax_adj_simple
- `extract_PR_constraints.m` - import results of optimisation at Cc = 129 umol mol-1, average the results and export PR enzyme Vmax values to serve as constraints for full optimizations
- `extract_optimization_results.m` - read in optimisation results for each Cc (130-380 umol mol-1) and average the results
- `hex2rgb.m` - convert hexadecimal colour code to RGB
- `import_optimization_results.m` - imports results data for non-optimised and optimised enzyme concentrations 
- `plot_A_Cc_new.m` - plot Cc vs A for FvCB and e-Photosynthesis models on a single graph
- `plot_calibration.m` - plot and examine results of the parameter calibration to identify the best scaling factors αEnzymes and αRubio
- `plot_calibration_Jmax.m` - plot and examine results of the parameter calibration to identify optimal αEnzymes
- `plot_calibration_Vcmax.m` - plot and examine results of the parameter calibration to identify optimal αRubisco
- `plot_models_2.m` - calculate assimilation for FvCB and e-Photosynthesis models
- `plot_scaling_factors.m` - plot sum-of-squared-residual (SSR) values and scaling factors obtained from Jmax_simple and Vcmax_simple
- `plot_violins_new.m` - create violin plots for assimilation strategies
- `Violin.m` - class file for violin plots
- `plotting_A_rates.m` - plot assimilation rates for each Cc (130-380) for non-optimised and optimised enzyme concentrations as well as applying assimilation strategies (stress, current, future)
- `plotting_bars_in_log_new_2.m` - plot bar charts for average non-optimised vs optimised enzyme concentrations in a loop
- `plotting_ind_scatter_averages.m` - create scatter plots for average non-optimised and optimised enzyme concentrations as individual enzyme plots
- `plotting_models.m` - plot scaling factors for Rubisco and RuBP-regeneration limitation using results from plot_models_2
- `plotting_scatter_averages_new.m` - create scatter plots for average non-optimised and optimised enzyme concentrations in a loop, creating panels for individual enzymes
- `plotting_scatter_reps_in_log.m` - create scatter plots for average optimised enzyme concentrations in a loop, according to functional categories of enzymes (CBB,PR,SS)
- `plotting_semilogy_all.m` - create semilogy plots for average non-optimised and optimised enzyme concentrations in a loop, creating panels for individual enzymes
- `plotting_semilogy_ind.m` - create semilogy plots for average non-optimised and optimised enzyme concentrations in a loop, according to functional categories of 7-8 enzymes (CBB,PR,SS)

Shell_Scripts - batch job scripts (.sh) for running `gpmain` simulations on the High-End Computing (HEC) cluster: https://lancaster-hec.readthedocs.io/
- simulations can be run by submitting single scripts or arrays (to run the same program multiple times)
- running `job_array_gpmain_rice_129.sh` creates the photorespiratory constraints for all other simulations

gpmain - scripts called by batch jobs for running model optimisation simulations at Cc = 130-380 umol mol-1
