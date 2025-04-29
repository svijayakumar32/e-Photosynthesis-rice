### e-Photosynthesis-rice
This repository contains the code and data necessary to reproduce the results presented in the paper by Vijayakumar, S., Wang, Y.,Lin, H.C., Carmo-Silva, E., Long, S.P. and Taylor, S. H. entitled "Tailoring a dynamic model of photosynthetic metabolism towards greater carbon assimilation in rice".

## Overview
This paper presents a computational workflow for optimising resource investment among photosynthetic enzymes to maximise carbon assimilation in a model of C3 photosynthetic metabolism calibrated to rice, summarised in `pipeline_simple_flowchart.pdf`.
This new model of rice integrates species-specific temperature dependencies for the balance of oxygenation and carboxylation (Vomax/Vcmax) and ribulose-1,6-bisphosphate carboxylase/oxygenase (Rubisco) activation, as well as accounting for tight binding Rubisco inhibitors.
The workflow begins by combining species-specific equations for calculating temperature dependencies of Rubisco catalytic properties with leaf-level gas exchange measurements for Oryza sativa cv. IR64 to derive photosynthetic parameters describing Calvin-Benson-Bassham (CBB) cycle activity, i.e. Vcmax and J. 
These photosynthetic parameters are used to re-scale enzyme activities in e-Photosynthesis before running the model to identify redistributions that are optimal for CO2 assimilation at different [CO2] levels. Following this, additional analyses are applied to evaluate strategies for increasing photosynthesis through overexpressing subsets of 2-6 enzymes to engineer improved photosynthesis under drought stress, current ambient [CO2] and future elevated [CO2].

## Requirements
Running the code requires installation of the following R packages:
- `devtools` - for installing other packages
- `msuRACiFit` - for fitting photosynthetic CO_2 response curves to assimilation 
- `here` - for constructing paths to project files
- `readxl` - to import Excel files into R

## Files
The files within the main repository contains scripts, variables and data associated with the version of e-Photosynthesis cloned from https://github.com/cropsinsilico/C3-metabolic-and-leaf-model.
- `cdn` defines the environmental conditions, such as CO2 concentration and photon flux density.
- `SYSInitial` defines the simulation time.
- `CM_Drive` describes carbon metabolism (CBB cycle).
- `PS_Drive` describes the CBB cycle in addition to the starch synthesis and triose phosphate export. 
- `PR_Drive` describes the photorespiration model. 
- `PS_PRDrive` includes the CBB cycle, starch synthesis, triose phosphate export, and photorespiration process. 
- `FI_Drive` includes the light energy absorption, transfer, primary charge separation, and electron transfer around PSII. 
- `BF_Drive` includes the electron transfer from reduced plastiquinone until the generation of ATP and NADPH, including the ion transfer through thylakoid membrane and ATP synthesis process. 
- `FIBF_Drive` includes all the reactions covered by FI_Drive, and BF_Drive. 
- `RuACT_Drive` includes reactions of Rubisco activation process. 
- `XanCycle_Drive` includes reaction of the xanthophylls cycle. 
- `EPS_Drive` includes all the reactions covered by PS_PRDrive and FIBF_Drive.
- `EPS_Drive_GRNs` calls EPS_Drive to calculate assimilation rates using a specific list of enzyme activities as input.
- `RA_Drive` includes reactions covered by EPS_Drive and RA_Drive. 
- `RedoxReg_Drive` includes reactions covered by EPS_Drive and the redox regulation of enzyme activities. 
- `DynaPS_Drive` includes reactions covered by EPS_Drive, RuACT_Drive and XanCycle_Drive.
- `tr DynaPS_Drive` includes reactions covered by EPS_Drive, RuACT_Drive, XanCycle_Drive and RROEA_Drive.

In certain filenames, the names of the drives are extended with suffixes:
- _Ini defines initial values and parameters.
- _Rate defines rate equations.
- _mb defines differential equations.

Additional custom scripts and data files are described below:
- `Calculate_Cc.m` - calculates Cc values equivalent to Ca for the FvCB model fit
- `CalculateGrossAssimilation.m` - calculates gross CO2 assimilation (Net_A + R_d)
- `CalculateNetAssimilation.m` - calculates net CO2 assimilation (Gross_A - R_d)
- `Einput7.txt` - a list of the original enzyme activities for e-Photosynthesis model

The sub-folders listed below contain files associated with each analysis:
# `CO_2 response` - 
# `Live_Scripts` - 
# `Outputs` - 
# `Parameterisation` - 
# `Plotting` - 
# `Shell_Scripts` -
# `gpmain` - 
