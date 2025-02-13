% Import data from excel file containing results of e-Photosynthesis optimization using rice data
file = 'Results_optimization_rice_new.xlsx';

% Define range to import data from
range = 'V5:X30'; %Average, STDEV, SEM for optimized enzyme protein concentrations

% Define sheet names to import data from
sheet_names = {'140','160','180','200','220','240','260','280','300','320','340','360','380','400','420'};

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

% Get enzyme names 
cats = readcell(file, 'Sheet', '140', 'Range', 'C5:C30'); % Read enzyme names
cats = categorical(cats(rows_to_keep, :)); % Remove duplicate enzymes and convert into categorical 

% Assign x and y for bar graphs
y = 1:23; % Y-axis values refer to enzymes
x = horzcat(results_data{1,1},results_data{2,15}(:,1)); % X-axis values refer to protein levels, replace second column results_data {2,i} for each Ci or make a loop e.g. 
error_values = results_data{2,15}(:,3); % Error values, replace for each Ci

% Convert data and error values to log scale
x_log = log10(x); % Log-transformed data
%error_values_log_old = log10(x + error_values) - x_log; % Log-transformed errors for non-opt and opt
error_values_log = log10(x(:,2) + error_values) - x_log(:,2); % Log-transformed errors for opt only

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

% ERROR BAR PLACEMENT - WORK IN PROGRESS
% % Calculate bar widths based on the difference between adjacent x-values
% bar_widths = abs(diff(x_log));
% 
% % Add error bars aligned with the x-endpoints of the bars
% hold on;
% for j = 1:numel(y)
%     % Calculate the x-position for the error bars (right side)
%     x_position = x_log(j) + bar_widths(j)/2; % Centered on the bar
% 
%     % Calculate the scaled error values for log scale
%     scaled_error = error_values(j) / (10^x_log(j));
% 
%     % Plot error bars for optimized (second bar)
%     errorbar(x_position, y(j), scaled_error, 'horizontal', 'k', 'Marker', 'none');
% end

% Add legend
legend('Non-optimized', 'Optimized', 'Location', 'NorthOutside');
hold off;

xlabel('Relative protein content (mg m^{-2})');
ylabel('Enzymes');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 10)

% Print new image in pdf document
print('barh_420','-dpdf','-fillpage');
