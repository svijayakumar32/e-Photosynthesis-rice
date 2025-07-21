This folder provides an example of how outputs of the enzyme optimisation can be stored in the same repository for ease of import by the plotting scripts.
- `Results_optimization_rice_new.xlsx` - an Excel workbook containing all results of the model optimisations - non-optimised and optimised enzyme Vmax, their corresponding enzyme concentrations, and assimilation rates following optimisation.

There are also subfolders to save these outputs separately:
- `Enzymes` - BestMatrix gives the optimal distribution of Vmax values for 67 photosynthetic enzymes (V1-V59)
- `Graphs` - A location to export plots from the optimisation/sensitivity analysis
- `Metabolites` - dplot gives change in metabolite concentrations (which reach steady state at the end of the optimisation)
- `Workspaces` - MATLAB workspaces saved after running gpmain simulations
