% Import data from excel file containing results of e-Photosynthesis optimization using rice data - This version loops both the generation of data and the graphing process
file = 'Results_optimization_rice.xlsx';

% Define range to import data from
range = 'V5:X30'; %Average, STDEV, SEM for optimized enzyme protein concentrations

% Define sheet names to import data from
sheet_names=num2cell(140:20:420);

% Create cell array to import non-optimized and optimized protein data from each Ci/sheet
results_data = cell(2, numel(sheet_names));

% Loop through all sheets and import values into results_data
for i = 1:numel(sheet_names)
    % Read the specified range as a numeric array
    numeric_data_unoptimized = readmatrix(file, 'Sheet', sheet_names{i}, 'Range', 'AI5:AI30');
    %numeric_data_optimized = readmatrix(file, 'Sheet', '140', 'Range', 'AI5:AI30');
    numeric_data_optimized = readmatrix(file, 'Sheet', sheet_names{i}, 'Range', range);
    % Replace zero values with NaN
    numeric_data_unoptimized(numeric_data_optimized == 0) = NaN;
    numeric_data_optimized(numeric_data_optimized == 0) = NaN;
    % Remove rows containing all NaN values
    rows_to_keep = any(~isnan(numeric_data_optimized), 2);
    numeric_data_optimized = numeric_data_optimized(rows_to_keep, :);
    numeric_data_unoptimized = numeric_data_unoptimized(rows_to_keep, :);
    % Store data in cell array using sheet names as indices (1-15 = 140-420)
    results_data{1,i} = numeric_data_unoptimized;
    results_data{2,i} = numeric_data_optimized;
end

% Calculate percentage changes per Ci
x = zeros(numel(y),2); % X-axis values refer to protein levels
for j = 1:numel(y)
    i = 1:numel(sheet_names);
    x(i) = horzcat(results_data{1,1},results_data{2,i}(:,1)); %Replace second column in results_data for each Ci
    percentage_change(j,i) = (x(i,2)-x(i,1))/x(i,1)*100;
end