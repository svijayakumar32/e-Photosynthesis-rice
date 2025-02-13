%% Plot and examine results of the parameter calibration to identify the best scaling factors αEnzymes and αRubisco

load('plot_data_J.mat');
load('plot_data_V.mat');

% Create figure
figure;
hold on;

% Specify Cc and Gross Assimilation data for FvCB model, joining data from
% Vcmax and J calibration
FvCB_x = vertcat(Farq_Matrix_V(1:(size(numbers_V,1)), 5),Farq_Matrix_J((1:size(numbers_J,1)), 5));
FvCB_y = vertcat(Farq_Matrix_V(1:(size(numbers_V,1)), 4),Farq_Matrix_J((1:size(numbers_J,1)), 4));

% Plot FvCB points with markers and dashed lines
plot(FvCB_x, FvCB_y, 'Marker','square', ...   % '--o' specifies dashed line with circle markers
                     'MarkerFaceColor', 'k', ... % k=black markers and line
                     'MarkerEdgeColor', 'k', ...
                     'MarkerSize', 8, ...
                     'Color','k',...
                     'DisplayName', 'FvCB');

% For e_Photosynthesis data, we are plotting multiple lines, each associated with a scaling factor between 0.52-1.5 (n=50)
% Each scaling factor has 15 Cc values (data points) associated with Vcmax 
% and 35 Cc values (data points) associated with Jmax
chunkSize_V = size(numbers_V,1);
chunkSize_J = size(numbers_J,1);

% Generate different colors for each new scaling factor
%numChunks remains same for Vcmax and J so just define either of them
numChunks = ceil(size(ePhoto_Matrix_V, 1) / chunkSize_V);
colors = lines(numChunks);

% Convert scaling factors into string array for a data table/legend, including FvCB as first series
% Again, scaling labels are same for both sets of scaling factors, so pick from either V or J
% But Cc labels and the assimilation values must be joined from both datasets
labels = string(SSR_Matrix_V(:, 1));
all_labels = vertcat("FvCB", labels);
Cc_labels = string(vertcat(Farq_Matrix_V(1:chunkSize_V,5),Farq_Matrix_J(1:chunkSize_J,5)));
Cc_labels_numeric = vertcat(Farq_Matrix_V(1:chunkSize_V,5),Farq_Matrix_J(1:chunkSize_J,5));
ePhoto_table = array2table(vertcat(numbers_V,numbers_J), 'VariableNames', labels,RowNames = Cc_labels);

% Define the Cc ranges for splitting the lines
range1 = 1:chunkSize_V; % Vcmax limited points
range2 = chunkSize_V+1:size(Cc_labels_numeric, 1); % J limited points

% We are choosing to only plot a subset of scaling factors and then the two optimal values 
scaling_idx_dashes = 5:10:45; % 0.6, 0.8, 1.0, 1.2 1.4 as dashed lines for comparison
%scaling_idx_alpha = [23,24]; % 0.96, 0.98 as solid lines in red with diff marker shapes
scaling_idx_alpha = [22,23]; % 0.94, 0.96 as solid lines in red with diff marker shapes

scaling_idx = horzcat(scaling_idx_dashes,scaling_idx_alpha);

% Get a subset of labels
%short_labels = vertcat("FvCB",labels(scaling_idx));

% Plot Cc (x) against Gross Assimilation (y) for the selected scaling factors (dashed lines) and the optimal values (red solid lines)
for i = 1:length(scaling_idx)
    current_idx = scaling_idx(i);
    
    % Extract data for the current scaling factor
    chunk_y = table2array(ePhoto_table(:, current_idx));
    
    % Check if the current index is an alpha value (for solid red lines)
    if current_idx == 22
        % Plot solid red lines with plus markers for alpha = 0.94
        plot(Cc_labels_numeric(range1), chunk_y(range1), ...
            'Marker','+', ...
            'Color', 'r', ... 
            'MarkerFaceColor','r',...
            'MarkerEdgeColor','r',...
            'MarkerSize', 8, ...
            'DisplayName', 'α = 0.94'); 
        plot(Cc_labels_numeric(range2), chunk_y(range2), ...
            'Marker','+', ...
            'Color', 'r', ... 
            'MarkerFaceColor','r',...
            'MarkerEdgeColor','r',...
            'MarkerSize', 8, ...
            'HandleVisibility', 'off'); 
    elseif current_idx == 23
        % Plot solid red lines with cross markers for alpha = 0.96
        plot(Cc_labels_numeric(range1), chunk_y(range1), ...
            'Marker','x', ...
            'Color', 'r', ... 
            'MarkerFaceColor','r',...
            'MarkerEdgeColor','r',...
            'MarkerSize', 8, ...
            'DisplayName', 'α = 0.96'); 
        plot(Cc_labels_numeric(range2), chunk_y(range2), ...
            'Marker','x', ...
            'Color', 'r', ... 
            'MarkerFaceColor','r',...
            'MarkerEdgeColor','r',...
            'MarkerSize', 8, ...
            'HandleVisibility', 'off'); 
    else
        % Plot dashed lines with circle markers for all other scaling factors
        plot(Cc_labels_numeric(range1), chunk_y(range1), ...
            '--o', ...
            'Color', colors(i, :), ...
            'MarkerFaceColor',colors(i, :), ...
            'MarkerEdgeColor',colors(i, :), ...
            'MarkerSize', 8, ...
    'DisplayName', sprintf('α = %s', labels(current_idx)));  
        plot(Cc_labels_numeric(range2), chunk_y(range2), ...
            '--o', ...
            'Color', colors(i, :), ...
            'MarkerFaceColor',colors(i, :), ...
            'MarkerEdgeColor',colors(i, :), ...
            'MarkerSize', 8, ...
            'HandleVisibility', 'off');
    end
end

% Add vertical dashed line at midpoint between Rubisco and RuBP_regeneration limited Cc ranges
%midpoint = Cc_labels_numeric(15) + (Cc_labels_numeric(16) - Cc_labels_numeric(15)) / 2;
%xline(midpoint, '--k');

hold off;
% 
% % Set the x-axis limits from 0 to 700
xlim([0 700]);
%xticks([0:100:640]);
set(gca, 'FontSize', 20); 
% 
% % Add labels, title, and legend
xlabel('C_c (μmol mol^{-1})','FontSize',20);
ylabel('Gross Assimilation (μmol m^{-2} s^{-1})','FontSize',20');
%title('Scaling Factors for Rubisco and RuBP-Regeneration Limitation','FontSize',12');

% Add legend, including FvCB
legend('show', 'Location', 'eastoutside');

% % Print in pdf
set(gcf, 'PaperOrientation', 'landscape');
print(gcf,fullfile('Outputs/rice_params/graphs',"Scaling_Factors"),'-dpdf','-fillpage');
