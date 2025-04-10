% This script imports data for non-optimized vs optimized enzyme concentrations in a loop

% Import data from excel file containing results of e-Photosynthesis optimization using rice data 
file = 'Outputs/rice_params/Results_optimization_rice_new_2.xlsx';

% Define range to import data from
non_optimized_range = 'AI5:AI30'; %
optimized_range = 'V5:X30'; %import averages, STDEV and SEM

% Define sheet names to import data from
sheet_names = string(num2cell(130:10:380)); % if there are regular intervals between values

% Create cell array to import non-optimized and optimized protein data from each Cc on a different sheet
results_data = cell(2, numel(sheet_names));

% Loop through all sheets and import values into results_data
for i = 1:numel(sheet_names)
    % Read the specified range as a numeric array
    numeric_data_non_optimized = readmatrix(file, 'Sheet', sheet_names{i}, 'Range', non_optimized_range);
    %numeric_data_optimized = readmatrix(file, 'Sheet', '140', 'Range', 'AI5:AI30');
    numeric_data_optimized = readmatrix(file, 'Sheet', sheet_names{i}, 'Range', optimized_range);
    % Replace zero values with NaN
    numeric_data_non_optimized(numeric_data_optimized == 0) = NaN;
    numeric_data_optimized(numeric_data_optimized == 0) = NaN;
    % Remove rows containing all NaN values
    rows_to_keep = any(~isnan(numeric_data_optimized), 2);
    numeric_data_optimized = numeric_data_optimized(rows_to_keep, :);
    numeric_data_non_optimized = numeric_data_non_optimized(rows_to_keep, :);
    % Store data in cell array using sheet names as indices (1-26 = 130-380)
    results_data{1,i} = numeric_data_non_optimized;
    results_data{2,i} = numeric_data_optimized;
end

% Get enzyme names 
cats = readcell(file, 'Sheet', '130', 'Range', 'C5:C30'); % Read enzyme names
cats = categorical(cats(rows_to_keep, :)); % Remove duplicate enzymes and convert into categorical

% Get PR constraints from optimizations at 129 ppm
% PR_constraints_protein = PR_constraints;
% Get protein constraints
PR_protein_rows = readmatrix('Outputs/rice_params/Results_optimization_rice_new_2.xlsx','Sheet','129','Range','AA17:AJ23');
PR_constraints_protein = num2cell(mean(PR_protein_rows,2));
%writematrix(PR_constraints_protein,'PR_constraints_protein.txt');
