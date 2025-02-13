% This script loops scatter plots for average non-optimized and optimized enzyme
% concentrations in a loop, creating panels for individual enzymes

%Import all required data in a separate script
import_optimization_results; %import average of ten replicates
%import_optimization_reps; % import ten replicates

% Assign x and y for scatter graphs
% Assign output variables to non-optimized and optimized protein data in linear and log scale 
y = cell(1, numel(sheet_names));
y_log = cell(1, numel(sheet_names));
y_fold = cell(1, numel(sheet_names));
%y_log_fold = cell(1, numel(sheet_names)); redundant

for i = 1:numel(sheet_names)
    % Combine full data for non-optimized and optimized at each Cc
    y{1,i} = horzcat(results_data{1,1},results_data{2,i}(:,1));
    % Convert data to log scale
    y_log{1,i} = horzcat(log10(y{1,i}(:,1)),log10(y{1,i}(:,2))); % log x
    % Calculate y fold changes
    y_fold{1,i} = y{1,i}(:,2:end)./y{1,i}(:,1);
    % Calculate y_log fold changes - essentially the same as y fold since
    % the ratios of non-opt to opt are equal
    %y_log_fold{1,i} = y_log{1,i}(:,2:end)./y_log{1,i}(:,1);
end

% Convert sheet_names into x values 
% Ci = str2double(sheet_names);
Cc = str2double(sheet_names);
% Ca = Ci/0.7;
enzymes=string(cats);

% Output variables for all enzymes
enzyme_nopt_concs = zeros(size(cats,1),size(y,2)); 
all_opt_concs = zeros(size(cats,1),size(y,2)); 
all_opt_FC = zeros(size(cats,1),size(y,2)); 
all_opt_log_FC = zeros(size(cats,1),size(y,2)); 

CBB_indices = 1:8;
PR_indices = 10:16;
SS_indices = [9,17:23];

% Replace PR_constraints reaction IDs with enzyme names 
PR_blanks = cell(7,1);
PR_constraints_protein = cat(2,PR_blanks,PR_constraints_protein);
PR_constraints_protein(:,1) = cellstr(enzymes(PR_indices));

% Horizontally stack PR_limits values for each Cc by repeating col 2 of PR_constraints
PR_limits = cell2mat(repmat(PR_constraints_protein(:,2),1,26));% last no represents no. of Cc levels
%PR_limits = cell2mat(repmat(PR_constraints_protein(:,2),1,15));

% Loop through all sheets and extract the values for specific enzymes
for j = 1:length(y{1,i}(:,2)) % For all enzymes
    for i = 1:numel(Cc)
    enzyme_nopt_concs(j,i) = y{1,i}(j,1); % take the non-optimized data from col 1
    all_opt_concs(j,i) = y{1,i}(j,2); % take the optimized data from col 2
    all_opt_FC(j,i) = y_fold{1,i}(j,1); 
    %all_opt_log_FC(j,i) = y_log_fold{1,i}(j,1); 
    end
end

%Run fold_change avg here- but no need anymore since Einput_140 uses Rubisco and SBPase Vmax optimised at 280 ppm%

% Use a switch to choose whether to plot non-opt vs opt or fold changes for each enzyme panel
plottype = input('Enter type of plot (opt/FC): ', 's');

% Create CBB panel
    figure('Position', [100, 100, 300, 2000]); % [left, bottom, width, height];
    tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    % Plot scatter graphs of Cc (x) against absolute protein concentrations of individual enzymes 
    for j = 1:length(CBB_indices) % For each enzyme j 
        CBB_index = CBB_indices(j);
        nexttile
        switch (plottype)
        case 'opt'
        %semilogy(Cc, enzyme_nopt_concs(CBB_index,:), '.', 'MarkerSize', 15);
        scatter(Cc, enzyme_nopt_concs(CBB_index,:),'black','filled','o'); % Plot non-optimized series
        hold on
        %semilogy(Cc, all_opt_concs(CBB_index,:), '^', 'MarkerSize', 15);
        scatter(Cc, all_opt_concs(CBB_index,:),'green','filled','^'); % Plot optimized series
        hold off
        % Cc x axis limits
        xlim([min(Cc), max(Cc)]); % Cc x axis limits
        xlim ("padded");
        %xticks(Cc)
        xticks(130:50:380);
        set(gca, 'FontSize', 10); 
        %xtickangle(90);
        %xticklabels(Cc);
        xticklabels(130:50:380);
        xlabel('C_c (μmol mol^{−1})','FontSize', 8);
        % Protein y axis limits 
        % % Automatically determine y ticks
        yticks_auto = get(gca, 'YTick');
        %yticks([0 yticks_auto]);
        %yticklabels(sprintf('%d\n', [0, yticks_auto]));
        % Label y axis
        ylabel('Average absolute protein content (mg m^{-2})', 'FontSize', 8);
        ylim padded
        % % Calculate tick interval
        tick_interval = mode(diff(yticks_auto));
        % % Determine max default tick intervals
        max_tick_value = max(yticks_auto);
        % Set upper limit of y-axis to be one interval higher than max
        ylim([0, max_tick_value + tick_interval]);
        %yticks([0, max_tick_value + tick_interval]);
        %yticklabels(sprintf('%d\n', [0, max_tick_value + tick_interval]));
        % Other properties
        title(enzymes(CBB_index))
        grid on
        hold off
            if j==length(CBB_indices)
            legend('Non-optimized','Optimized','Location', 'SouthEastOutside')
            end
    %%%OR only one series for fold change
        case 'FC'
        scatter(Cc, all_opt_FC(CBB_index,:),'blue','filled','^'); % Plot fold change series
        hold on
        xlim([min(Cc), max(Cc)]); % Cc x axis limits
        ylim([0,inf]) % Protein y axis limits
        yticklabels(0, inf)
        %xticks(Cc)
        xticks(130:50:380);
        %xticklabels(Cc)
        xticklabels(130:50:380);
        xlabel('C_c (μmol mol^{−1})','FontSize', 10);
        ylabel('Fold change protein content','FontSize', 10);
        title(enzymes(CBB_index))
        grid on
        hold off
        end
    end
set(gcf, 'PaperOrientation', 'portrait');
print(gcf,fullfile('Outputs/rice_params/graphs',"Cc_vs_CBB_Enzymes"),'-dpdf','-fillpage');

% Create PR panel
figure('Position', [100, 100, 300, 2000]); % [left, bottom, width, height];
tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% Plot scatter graphs of Cc (x) against absolute protein concentrations of individual enzymes 
for j = 1:length(PR_indices) % For each enzyme j 
    PR_index = PR_indices(j);
    nexttile
    switch (plottype)
        case 'opt'
        %semilogy(Cc, enzyme_nopt_concs(PR_index,:), '.', 'MarkerSize', 15);
        scatter(Cc, enzyme_nopt_concs(PR_index,:),'black','filled','o'); % Plot non-optimized series
        hold on
        %semilogy(Cc, all_opt_concs(PR_index,:), '^', 'MarkerSize', 15);
        scatter(Cc, all_opt_concs(PR_index,:),'green','filled','^'); % Plot optimized series
        hold on 
        % Plot PR constraints levels
        plot(Cc, PR_limits(j,:), "b--"); 
        hold off
        % Cc x axis limits
        xlim([min(Cc), max(Cc)]); % Cc x axis limits
        xlim ("padded");
        %xticks(Cc)
        xticks(130:50:380);
        set(gca, 'FontSize', 10); 
        %xtickangle(90);
        %xticklabels(Cc)
        xticklabels(130:50:380);
        xlabel('C_c (μmol mol^{−1})','FontSize', 8);
        % Protein y axis limits 
        % % Automatically determine y ticks
        yticks_auto = get(gca, 'YTick');
        %yticks([0 yticks_auto]);
        %yticklabels(sprintf('%d\n', [0, yticks_auto]));
        % Label y axis
        ylabel('Average absolute protein content (mg m^{-2})', 'FontSize', 8);
        ylim padded
        % % Calculate tick interval
        tick_interval = mode(diff(yticks_auto));
        % % Determine max default tick intervals
        max_tick_value = max(yticks_auto);
        % Set upper limit of y-axis to be one interval higher than max
        ylim([0, max_tick_value + tick_interval]);
        title(enzymes(PR_index))
        grid on
        hold off
            if j==length(PR_indices)
            legend('Non-optimized','Optimized','Lower Limit','Location','SouthEastOutside')
            end
        %%%OR only one series for fold change
        case 'FC'
        scatter(Cc, all_opt_FC(PR_index,:),'blue','filled','^'); % Plot fold change series
        hold on
        xlim([min(Cc), max(Cc)]); % Cc x axis limits
        ylim([0,inf]) % Protein y axis limits 
        %xticks(Cc)
        xticks(130:50:380);
        %xticklabels(Cc)
        xticklabels(130:50:380);
        xlabel('C_c (μmol mol^{−1})','FontSize', 10);
        ylabel('Fold change protein content','FontSize', 10);
        title(enzymes(PR_index))
        grid on
        hold off
     end
end
set(gcf, 'PaperOrientation', 'portrait');
print(gcf,fullfile('Outputs/rice_params/graphs',"Cc_vs_PR_Enzymes"),'-dpdf','-fillpage');

% Create SS panel
figure('Position', [100, 100, 300, 2000]); % [left, bottom, width, height];
tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% Plot scatter graphs of Cc (x) against absolute protein concentrations of individual enzymes 
for j = 1:length(SS_indices) % For each enzyme j 
    SS_index = SS_indices(j);
    nexttile
        switch (plottype)
        case 'opt'
        %semilogy(Cc, enzyme_nopt_concs(SS_index,:), '.', 'MarkerSize', 15);
        scatter(Cc, enzyme_nopt_concs(SS_index,:),'black','filled','o'); % Plot non-optimized series
        hold on
        %semilogy(Cc, all_opt_concs(SS_index,:), '^', 'MarkerSize', 15);
        scatter(Cc, all_opt_concs(SS_index,:),'green','filled','^'); % Plot optimized series
        hold off
        % Cc x axis limits
        xlim([min(Cc), max(Cc)]); % Cc x axis limits
        xlim ("padded");
        %xticks(Cc)
        xticks(130:50:380);
        set(gca, 'FontSize', 10); 
        %xtickangle(90);
        %xticklabels(Cc)
        xticklabels(130:50:380);
        xlabel('C_c (μmol mol^{−1})','FontSize', 8);
        % Protein y axis limits 
        % % Automatically determine y ticks
        yticks_auto = get(gca, 'YTick');
        %yticks([0 yticks_auto]);
        %yticklabels(sprintf('%d\n', [0, yticks_auto]));
        % Label y axis
        ylabel('Average absolute protein content (mg m^{-2})', 'FontSize', 8);
        ylim padded
        % % Calculate tick interval
        tick_interval = mode(diff(yticks_auto));
        % % Determine max default tick intervals
        max_tick_value = max(yticks_auto);
        % Set upper limit of y-axis to be one interval higher than max
        ylim([0, max_tick_value + tick_interval]);
        title(enzymes(SS_index))
        grid on
        hold off
            if j==length(SS_indices)
            legend('Non-optimized','Optimized','Location', 'SouthEastOutside')
            end        %%%OR only one series for fold change
        case 'FC'
        scatter(Cc, all_opt_FC(SS_index,:),'blue','filled','^'); % Plot fold change series
        hold on
        xlim([min(Cc), max(Cc)]); % Cc x axis limits
        ylim([0,inf]) % Protein y axis limits 
        %xticks(Cc)
        xticks(130:50:380);
        %xticklabels(Cc)
        xticklabels(130:50:380);
        xlabel('C_c (μmol mol^{−1})','FontSize', 10);
        ylabel('Fold change protein content','FontSize', 10);
        title(enzymes(SS_index))
        grid on
        hold off
        end
end
set(gcf, 'PaperOrientation', 'portrait');
print(gcf,fullfile('Outputs/rice_params/graphs',"Cc_vs_SS_Enzymes"),'-dpdf','-fillpage');