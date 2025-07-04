### e-Photosynthesis-rice
This repository contains the code and data necessary to reproduce the results presented in the paper by Vijayakumar, S., Wang, Y.,Lin, H.C., Carmo-Silva, E., Long, S.P. and Taylor, S. H. entitled "Tailoring a dynamic model of photosynthetic metabolism towards greater carbon assimilation in rice".

## Overview
This paper presents a computational workflow for optimising resource investment among photosynthetic enzymes to maximise carbon assimilation in a model of C3 photosynthetic metabolism calibrated to rice, summarised in `pipeline_simple_flowchart.png`.
This version of e-Photosynthesis integrates species-specific temperature dependencies for the balance of oxygenation and carboxylation (Vomax/Vcmax) and ribulose-1,6-bisphosphate carboxylase/oxygenase (Rubisco) activation, as well as accounting for tight binding Rubisco inhibitors.
The workflow begins by combining species-specific equations for calculating temperature dependencies of Rubisco catalytic properties with leaf-level gas exchange measurements for Oryza sativa cv. IR64 to derive photosynthetic parameters describing Calvin-Benson-Bassham (CBB) cycle activity, i.e. Vcmax and J. 
These photosynthetic parameters are used to re-scale enzyme activities in e-Photosynthesis before running the model to identify redistributions that are optimal for CO2 assimilation at different [CO2] levels. Following this, additional analyses are applied to evaluate strategies for increasing photosynthesis through overexpressing subsets of 2-6 enzymes to engineer improved photosynthesis under drought stress, current ambient [CO2] and future elevated [CO2].

## Requirements
Running the code requires installation of the following R packages:
- `devtools` - for installing other packages.
- `msuRACiFit` - for fitting photosynthetic CO2 response curves to assimilation. 
- `here` - for constructing paths to project files.
- `readxl` - to import Excel files into R.

## Files
The files within the main repository contains scripts, variables and data adapted from the version of e-Photosynthesis stored in https://github.com/cropsinsilico/C3-metabolic-and-leaf-model.
All code licensing information is detailed in `LICENSE`.

- `cdn.m` defines the environmental conditions, such as CO2 concentration and photon flux density.
- `SYSInitial.m` defines the simulation time.

The evolutionary algorithm is run using the following scripts:
- `average.m` - calculates the average Vmax in each generation
- `mutate.m` - generate variation in Vmax populations by modifying mutatePercentage 
- `rankPop.m` - ranks populations based on their CO2 uptake values
- `resizePop.m` - resizes the population after each generation
- `stdev.m` - calculate the standard deviation of Vmax in each generation
- `twoPoint.m` - a variation of `resizePop` (not used in current simulations)
  
For the model scripts, the first part of the filename describes the portion of the e-Photosynthesis model being simulated:
- `CM_` describes carbon metabolism (CBB cycle).
- `PS_` describes the CBB cycle in addition to the starch synthesis and triose phosphate export. 
- `PR_` describes the photorespiration model. 
- `PS_PR_` describes the CBB cycle, starch synthesis, triose phosphate export, and photorespiration process. 
- `FI_` describes the light energy absorption, transfer, primary charge separation, and electron transfer around PSII. 
- `BF_` describes the electron transfer from reduced plastiquinone until the generation of ATP and NADPH, including the ion transfer through thylakoid membrane and ATP synthesis process. 
- `FIBF_` describes reactions covered by FI_Drive, and BF_Drive. 
- `RuACT_` describes reactions of Rubisco activation process. 
- `XanCycle_` describes reaction of the xanthophyll cycle. 
- `EPS_` describes reactions covered by PS_PRDrive and FIBF_Drive.
- `RA_` describes reactions covered by EPS_Drive and RA_Drive. 
- `RedoxReg_` describes reactions covered by EPS_Drive and the redox regulation of enzyme activities. 
- `DynaPS_` describes reactions covered by EPS_Drive, RuACT_Drive and XanCycle_Drive.
- `tr DynaPS_` describes reactions covered by EPS_Drive, RuACT_Drive, XanCycle_Drive and RROEA_Drive.

The second part of these filenames use suffixes to describe the purpose of the script:
- `_AddTitle` defines titles for concentration/time plots plots.
- `_Ini` defines initial values and parameters.
- `_Drive` defines model simulation settings.
- `_Graph` defines plotting commands.
- `_Rate` defines rate equations.
- `_mb` defines differential equations.

Custom scripts and data files are described below:
- `Calculate_Cc.m` - calculates Cc values equivalent to Ca for the FvCB model fit
- `CalculateGrossAssimilation.m` - calculates gross CO2 assimilation
- `CalculateNetAssimilation.m` - calculates net CO2 assimilation
- `Einput7.txt` - a list of the original enzyme activities for e-Photosynthesis model
- `EPS_Drive_GRNs` - calls EPS_Drive to calculate assimilation rates using a specific list of enzyme activities as input.
- `enzyme_adjustment_test_new`- runs sensitivity analysis to determine assimilation rates associated for random fold changes (n = 2000) of various enzyme combinations drawn between FC = 1 (no change) and FC = 1.25 for Rubisco and FC = 2.8 for SBPase, aldolase, PRK, FBPase and TK at Cc = 130 umol mol-1, 250 umol mol-1 or Cc = 360 umol mol-1
- `FvCB_check_Vcmax_adj.m` - checks transition between Rubisco-limited and RuBP limited phases in A/Cc fit
- `IModelCom.m` - called by CM_Drive to initialise the structure of the model using different components of the full photosynthesis model
- `IniModelCom.m` - called by EPS_Drive to to initialise the structure of the model using different components of the full photosynthesis model
- `Jmax_adj_simple.m` - runs model calibration to find optimal αEnzymes by minimising SSR between assimilation rates of FvCB and e-Photosynthesis models in the RuBP-regeneration limited range of [CO2]
- `MW&Kcat.txt` - a list of molecular weights and kcat values for photosynthetic enzymes
- `PR_constraints.txt` - a list of photorespiratory enzyme activities (optimised Vmax) obtained from running `job_array_gpmain_rice_129.sh`
- `PR_constraints_protein.txt` - a list of photorespiratory enzyme concentrations (optimised protein) converted from `PR_constraints`
- `Vcmax_adj_simple.m` - runs model calibration to find optimal αRubisco by minimising SSR between assimilation rates of FvCB and e-Photosynthesis models in the Rubisco limited range of [CO2]
  
In addition to this, the remaining files are stored in various directories defining their function:

CO_2 response - for fitting leaf-level gas exchange measurements to the Farquhar-vonCaemmerer-Berry (FvCB) model
- `fitting_gas_exchange_rice.Rproj` - R project

The Data sub-folder contains:
- `Gas_exchange_measurement_WT_plants.xls` - Excel workbook containing 8 sets of leaf level gas exchange measurements for wild type Oryza sativa cv. IR64
- Filenames beginning with `IR64-A009-07-33-05` contain data for each curve 
 
The Scripts sub-folder contains:
- `msuRACiFit_rice` - R Markdown document containing fitting procedure for obtaining photosynthetic parameters

Live_Scripts - for running interactive executable versions of scripts - combining code, output, and formatted text in a single MATLAB environment

Outputs - for storing outputs from various analyses within the subfolders:
- `Data` - contains `Optimization_Procedure.txt`, `Results_optimization_rice_new_2.xlsx`, and `Series_of_job_submissions.txt`
- `Enzymes` - BestMatrix gives the optimal distribution of Vmax values for 67 photosynthetic enzymes (V1-V59)
- `Metabolites` - dplot gives change in metabolite concentrations (which reach steady state at the end of the optimisation)
- `Workspaces` - MATLAB workspaces saved after running gpmain simulation

Parameterisation - additional functions for calculating model parameters
- `calculate_O2_sol` - to calculate molar solubility of oxygen at a given temperature and pressure
- `calculate_PsPr_ratio` - to calculate the Vo/Vc ratio for Rubisco based on temperature
- `find_Rubisco_params` - to find the value of the Rubisco rate parameter PsV1 at different temperatures

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
