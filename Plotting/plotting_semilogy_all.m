% This script loops scatter plots for average optimized enzyme concentrations in a
% loop, according to functional categories of 7-8 enzymes (CBB,PR,SS)

%Import all required data in a separate script
import_optimization_results; %import average of ten replicates
%import_optimization_reps; % import ten replicates

% Assign x and y for scatter graphs
% Assign output variables to non-optimized and optimized protein data in linear and log scale 
y = cell(1, numel(sheet_names));
y_log = cell(1, numel(sheet_names));
y_fold = cell(1, numel(sheet_names));

for i = 1:numel(sheet_names)
    % Combine full data for non-optimized and optimized at each Ci
    y{1,i} = horzcat(results_data{1,1},results_data{2,i}(:,1));
    % Convert data to log scale
    y_log{1,i} = horzcat(log10(y{1,i}(:,1)),log10(y{1,i}(:,2))); % log x
end

% Convert sheet_names into x values 
Ci = str2double(sheet_names);
Ca = Ci/0.7;

% Specify enzyme of interest
%rubisco_concs = zeros(size(y_log));
%enzyme_nopt_concs = zeros(size(cats,1),size(y,2)); 

% This is kinda pointless for the  linear values since we can just subset results_data but we need y_log values for log10 of these values
CBB_opt_concs = zeros(8,size(y,2)); % CBB cycle enzymes
PR_opt_concs = zeros(7,size(y,2));  % Photorespiratory enzymes
SS_opt_concs = zeros(8,size(y,2));  % Sucrose/starch enzymes

CBB_FC = zeros(8,size(y,2)); % CBB cycle enzymes
PR_FC = zeros(7,size(y,2));  % Photorespiratory enzymes
SS_FC = zeros(8,size(y,2));  % Sucrose/starch enzymes

CBB_indices = 1:8;
PR_indices = 10:16;
SS_indices = [9,17:23];

% Loop through all sheets and extract the values for specific enzyme groups

% e.g. All CBB enzymes
for j = 1:length(CBB_indices) % Rubisco - PRK
    CBB_index = CBB_indices(j);
    for i = 1:numel(Ci)
    %enzyme_nopt_concs(j,i) = y{1,i}(j,1); 
    CBB_opt_concs(j,i) = y{1,i}(CBB_index,2); 
    CBB_FC(j,i) = y{1,i}(CBB_index,2)./y{1,i}(CBB_index,1); 
    end
end

% e.g. All PR enzymes
for j = 1:length(PR_indices) % PGP - GDC
    PR_index = PR_indices(j);
    for i = 1:numel(Ci)
    %enzyme_nopt_concs(j,i) = y{1,i}(j,1); 
    PR_opt_concs(j,i) = y{1,i}(PR_index,2); 
    PR_FC(j,i) = y{1,i}(PR_index,2)./y{1,i}(PR_index,1);
    end
end

% e.g. All SS enzymes
for j = 1:length(SS_indices) % AGPase, Aldolase(C) - 6P2K
    SS_index = SS_indices(j);
    for i = 1:numel(Ci)
    %enzyme_nopt_concs(j,i) = y{1,i}(j,1); 
    SS_opt_concs(j,i) = y{1,i}(SS_index,2); 
    SS_FC(j,i) = y{1,i}(SS_index,2)./y{1,i}(SS_index,1);
    end
end

%Specify a new colormap
%'forestgreen','plum','royalblue','burntorange','orangered','limegreen','spirodiscoball'

% Define hexadecimal color codes
%hex_colors = {'#014421','#8E4585','#002366','#CC5500','#FF4500','#32CD32','0FC0FC'};

figure;
grid on
nexttile
% Plot scatter graph of Ci (x) and log protein concentration of enzyme e.g. Rubisco
for j = 1:length(CBB_indices) % Rubisco - PRK - each j should be a new series/enzyme
    semilogy(Ci, CBB_opt_concs(j,:), '.', 'MarkerSize', 20);
    xlim([min(Ci), max(Ci)]); % Ci axis limits
    xticks(Ci)
    xticklabels(Ci)
    xlabel('C_i (ppm)');
    ylabel('Average absolute protein content (mg m^{-2})');
    title('Calvin-Benson-Bassham Cycle Enzymes')
    hold on
    % Legend
    legend(cats(CBB_indices),'Location', 'EastOutside')
end

% figure;
nexttile
grid on
% Plot scatter graph of Ci (x) and log protein concentration of enzyme 
for j = 1:length(PR_indices) 
    semilogy(Ci, PR_opt_concs(j,:), '.', 'MarkerSize', 20);
    xlim([min(Ci), max(Ci)]); % Ci axis limits
    xticks(Ci)
    xticklabels(Ci)
    xlabel('C_i (ppm)');
    ylabel('Average absolute protein content (mg m^{-2})');
    title('Photorespiratory Enzymes')
    hold on
    % Legend
    legend(cats(PR_indices),'Location', 'EastOutside')
end

% figure;
nexttile
grid on
% Plot scatter graph of Ci (x) and absolute protein concentration of enzyme 
for j = 1:length(SS_indices) 
    semilogy(Ci, SS_opt_concs(j,:), '.', 'MarkerSize', 20);
    xlim([min(Ci), max(Ci)]); % Ci axis limits
    xticks(Ci)
    xticklabels(Ci)
    xlabel('C_i (ppm)');
    ylabel('Average absolute protein content (mg m^{-2})');
    title('Sucrose/Starch Enzymes')
    hold on
    % Legend
    legend(cats(SS_indices),'Location', 'EastOutside')
end

print(fullfile('Outputs/rice_params/graphs','Ci_vs_Opt_Concs'),'-dpdf','-fillpage');
