%% Plotting bar charts for non-optimized vs optimized enzyme concentrations
% This function returns a bar chart taking inputs of the filename,
% sheet_name i.e. Ci and the ranges (cell_references) of non-optimized and optimized concentration values
function [results_data,x,x_log,error_values,error_values_log,h] = plot_bars_in_log(file,sheet_name,non_optimized_range,optimized_range)

% Create cell array to import non-optimized and optimized protein data from each Ci/sheet
results_data = cell(2, 1);

% Import values into results_data
% Read the enzyme names as a cell array
cats = readcell(file, 'Sheet', sheet_name, 'Range', 'C5:C30'); % Read enzyme names
% Read the specified range as a numeric array
numeric_data_non_optimized = readmatrix(file, 'Sheet', sheet_name, 'Range', non_optimized_range);
numeric_data_optimized = readmatrix(file, 'Sheet', sheet_name, 'Range', optimized_range);
% Replace zero values with NaN
numeric_data_non_optimized(numeric_data_optimized == 0) = NaN;
numeric_data_optimized(numeric_data_optimized == 0) = NaN;
% Remove rows containing all NaN values
rows_to_keep = any(~isnan(numeric_data_optimized), 2);
numeric_data_optimized = numeric_data_optimized(rows_to_keep, :);
numeric_data_non_optimized = numeric_data_non_optimized(rows_to_keep, :);
% Remove enzyme rows containing all NaN values and convert into categorical 
cats = categorical(cats(rows_to_keep, :)); 
% Store data in cell array using sheet names as indices (1-15 = 140-420)
results_data{1,1} = numeric_data_non_optimized;
results_data{2,1} = numeric_data_optimized;

% Assign x and y for bar graphs
y = 1:23; % Y-axis values refer to enzymes, doesn't change
% Assign output variables to non-optimized and optimized protein data  
x = horzcat(results_data{1,1},results_data{2,1}(:,1)); %Concatenate unoptimized and optimized protein levels
error_values = results_data{2,1}(:,3); % Error values
% Convert data and error values to log scale
x_log = log10(x); % Log-transformed protein data
error_values_log = log10(x(:,2) + error_values) - x_log(:,2); % Log-transformed SEMs of X for optimized only

% Plot data of optimized vs non-optimized enzyme concentrations
figure;
h = barh(y, x_log);

% Specify appearance of horizontal bars
h(1).FaceColor = "#000000"; % Black for non-optimized
h(2).FaceColor = "#92D050"; % Green for optimized
h(1).LineWidth = 0.25; % Reduce thickness of outlines
h(2).LineWidth = 0.25;

% Plot categories in reverse, manually adjust values in log scale, and use labels with linear values
set(gca, 'YDir', 'reverse', 'YTick', 1:23, 'YTickLabel', cats, 'Xtick', -2:4, 'Xticklabel', 10.^get(gca, 'Xtick'))

% Add legend
legend('Non-optimized', 'Optimized', 'Location', 'NorthOutside');
hold off;

xlabel('Relative protein content (mg m^{-2})');
ylabel('Enzymes');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 10)

% Print new image in pdf document
%graph_filename = strcat ("barh_",sheet_name);    
%print(fullfile('Outputs/rice_params/graphs',graph_filename),'-dpdf','-fillpage');

end