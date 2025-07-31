## Structure

This repository contains code for fitting CO₂ response data using the `msuRACiFit` package.

The analysis is fully scripted within a single .R file for portability and reproducibility.

All core logic is contained in `msuRACiFit_rice_params.R`, located in the `Scripts/` folder.

The `Data/' subdirectory contains the gas exchange data files used in the analysis.

- Raw data:
`Gas_exchange_measurement_WT_plants.csv`

- Cleaned input files:
These are named using the pattern:
`IR64-A009-07-33-05-0x_Wildtypex.csv`

## Requirements
Running the code requires installation of the following R packages:
- `devtools` - for installing other packages.
- `msuRACiFit` - for fitting photosynthetic CO2 response curves to assimilation. 
- `here` - for constructing paths to project files.
- `readxl` - to import Excel files into R.

## Running the Analysis

To run the full analysis:

1. Open a terminal or command prompt.

2. Navigate to the repository root directory.

3. Run the script using one of the following approaches:

     (a) Modify the setwd() line in `msuRACiFit_rice_params.R` to match your local path to the repository, then run:
     `source("Scripts/x_params.R")`
   
     (b) Alternatively, install the here package using:
     `install.packages("here")`.
   
   Then replace the `setwd()` line with:
   `library(here)`
   and run the script using:
   `source(here("Scripts", "msuRACiFit_rice_params.R"))`.

This allows the script to run from any environment within the directory.
