% This script loops scatter plots for average optimized enzyme concentrations in a
% loop, according to functional categories of 7-8 enzymes (CBB,PR,SS)

%Import all required data in a separate script
import_optimization_reps; % import ten replicates

% Convert sheet_names into x values 
Ci = str2double(sheet_names);
Ca = Ci/0.7;

% Assign output variables to non-optimized and optimized protein data in linear and log scale  for y values
y = cell(1, numel(sheet_names));
y_log = cell(1, numel(sheet_names));
y_fold = cell(1, numel(sheet_names));
y_log_fold = cell(1, numel(sheet_names));
mean_data = cell(length(results_data{1,1}), numel(sheet_names));
sem_data = cell(length(results_data{1,1}), numel(sheet_names));

for i = 1:numel(sheet_names) %15
    for j =1:length(results_data{1,1}) %23
    % Combine full data for non-optimized and optimized at each rep and Ci
    y{1,i} = horzcat(results_data{1,1},results_data{2,i}(:,1:end));
    % Convert data to log scale
    y_log{1,i} = horzcat(log10(y{1,i}(:,1)),log10(y{1,i}(:,2:end))); % log x
    % Calculate fold change between optimized/non-optimized at each Ci
    y_fold{1,i} = y{1,i}(:,2:end)./y{1,i}(:,1);
    % Calculate log fold change between optimized/non-optimized at each Ci
    y_log_fold{1,i} = y_log{1,i}(:,2:end)./y_log{1,i}(:,1);
    % Mean and standard error of means across ten replicates for each Ci
    mean_data{j,i} = mean(y{1,i}(j,2:10));
    sem_data{j,i} = std(y{1,i}(j,2:10))/sqrt(10);
    end
end

% Plot all reps for an enzyme with error bars
for j =1:length(y_fold{1,1})
    for i = 1:numel(sheet_names)
    errorbar(Ci(i), mean_data{j,i}, sem_data{j,i}, 'o-');  
    grid on;
    hold on;
    end
    hold off
    xlabel('C_i (ppm)');
    xticks(Ci)
    xticklabels(Ci)
    ylabel('Average absolute protein content (mg m^{-2})');
    enzymes=string(cats);
    title(enzymes(j));
    graph_filename = strcat ("Ci_vs_Reps_",enzymes(j));    
    print(fullfile('Outputs/rice_params/graphs/Ci_vs_Reps',graph_filename),'-dpdf','-fillpage');
end
