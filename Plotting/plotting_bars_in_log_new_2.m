% This script loops generation of bar charts for average non-optimized vs optimized enzyme concentrations in a loop, 
% returning all 15 graphs in series
% Alternatively, use the function plot_bars_in_log to generate graphs for each Ci at a time

%Import all required data in a separate script
import_optimization_results; 

% Assign x and y for bar graphs
y = 1:23; % Y-axis values refer to enzymes, doesn't change
% Assign output variable to non-optimized and optimized protein data  
x = zeros(numel(y),2); % X-axis values refer to protein levels
x_log = zeros(numel(y),2); % Logged protein levels
error_values = zeros(numel(y),1); % Corresponding SEMs of X 
error_values_log = zeros(numel(y),1); % Logged protein levels

for i = 1:numel(sheet_names)
    x = horzcat(results_data{1,1},results_data{2,i}(:,1)); %Replace second column in results_data for each Ci
    error_values = results_data{2,i}(:,3); % Error values, replace for each Ci
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

    % % ERROR BAR PLACEMENT - WORK IN PROGRESS
    % % Calculate bar widths based on the difference between adjacent x-values
    % bar_widths = abs(diff(x_log)); 
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
    graph_filename = strcat ("barh_",sheet_names(i));    
    print(fullfile('Outputs/rice_params/graphs',graph_filename),'-dpdf','-fillpage');
 end

% If calling the plot_bars_in_log function: type the following, replacing Ci in the sheetname 
% plot_bars_in_log('Results_optimization_rice.xlsx','140','AI5:AI30','V5:X30')
