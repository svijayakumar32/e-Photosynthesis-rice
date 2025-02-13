%%% Plot SSR and scaling factors obtained from Jmax_simple and Vcmax_simple
%% Load data
load Jmax_simple_new_result.mat;
load Vcmax_simple_new_result.mat;
%% Plot SSR against scaling factors and highlight the optimal points
% RuBP regeneration
figure;
plot(SSR_Matrix_J(:,1),SSR_Matrix_J(:,2),'g-o',...
     'MarkerFaceColor', '#5EBB50', ... % g = green markers and line
     'MarkerEdgeColor', '#5EBB50', ...
     'MarkerSize', 6, ...
     'Color','#5EBB50')
hold on
plot(SSR_Matrix_J(scaling_index_J,1),SSR_Matrix_J(scaling_index_J,2),...
            'o', ...
            'Color', 'r', ...
            'MarkerFaceColor', 'r', ... % r = red marker
            'MarkerEdgeColor', 'r', ...
            'MarkerSize', 8, ...
            'Color','r');
set(gca, 'xtick', 0.5:0.2:1.5);
set(gca, 'FontSize', 10);
xlabel('α_{Enzymes}','FontSize', 12)
ylabel('SSR','FontSize', 12)
hold off
% Export figure
%set(gcf, 'PaperOrientation', 'landscape');
%print(gcf,fullfile('Outputs/rice_params/graphs',"Scaling_Factors_for_RuBP_regeneration_new"),'-dpdf','-bestfit');
%% Rubisco carboxylation
figure;
plot(SSR_Matrix_V(:,1),SSR_Matrix_V(:,2),'g-o',...
     'MarkerFaceColor', '#5EBB50', ... % g = green markers and line
     'MarkerEdgeColor', '#5EBB50', ...
     'MarkerSize', 6, ...
     'Color','#5EBB50')
hold on
plot(SSR_Matrix_V(scaling_index_V,1),SSR_Matrix_V(scaling_index_V,2),...
            'o', ...
            'Color', 'r', ...
            'MarkerFaceColor', 'r', ... % r = red marker
            'MarkerEdgeColor', 'r', ...
            'MarkerSize', 8, ...
            'Color','r');
set(gca, 'xtick', 0.5:0.2:1.5);
set(gca, 'FontSize', 10);
xlabel('α_{Rubisco}','FontSize', 12)
ylabel('SSR','FontSize', 12)
hold off
% Export figure
%set(gcf, 'PaperOrientation', 'portrait');
%print(gcf,fullfile('Outputs/rice_params/graphs',"Scaling_Factors_for_Rubisco_carboxylation_new"),'-dpdf','-bestfit');

%% Create figure for RuBP regeneration scaling factors
figure;
% Specify Cc and Gross Assimilation data for FvCB model
FvCB_J_x = Farq_Matrix_J(1:4,5); % Cc - only first four rows needed as they repeat
FvCB_J_y = Farq_Matrix_J(1:4,4); % Gross A - only first four rows needed as they repeat

% Specify Cc and Gross Assimilation data for optimal scaling factor
% Round off scaling factor values to prevent errors associated with floating-point numbers
ePhoto_Matrix_J(:, 1) = round(ePhoto_Matrix_J(:, 1), 2);
% Find rows for scaling_factors in ePhoto by finding and assigning indices of ePhotoMatrix from Jmax_result

opt_indices_J = find(ePhoto_Matrix_J(:,1)==a_Enzymes);
idx06_J = find(ePhoto_Matrix_J(:,1)==0.6);
idx08_J = find(ePhoto_Matrix_J(:,1)==0.8);
idx10_J = find(ePhoto_Matrix_J(:,1)==1.0);
idx12_J = find(ePhoto_Matrix_J(:,1)==1.2);
idx14_J = find(ePhoto_Matrix_J(:,1)==1.4);

ePhoto_opt_x_J = ePhoto_Matrix_J(opt_indices_J,5); % Cc
ePhoto_opt_y_J = ePhoto_Matrix_J(opt_indices_J,4); % Gross A

ePhoto_J_x_06 = ePhoto_Matrix_J(idx06_J,5);
ePhoto_J_y_06 = ePhoto_Matrix_J(idx06_J,4);

ePhoto_J_x_08 = ePhoto_Matrix_J(idx08_J,5);
ePhoto_J_y_08 = ePhoto_Matrix_J(idx08_J,4);

ePhoto_J_x_10 = ePhoto_Matrix_J(idx10_J,5);
ePhoto_J_y_10 = ePhoto_Matrix_J(idx10_J,4);

ePhoto_J_x_12 = ePhoto_Matrix_J(idx12_J,5);
ePhoto_J_y_12 = ePhoto_Matrix_J(idx12_J,4);

ePhoto_J_x_14 = ePhoto_Matrix_J(idx14_J,5);
ePhoto_J_y_14 = ePhoto_Matrix_J(idx14_J,4);

% Plot FvCB points with markers and dashed lines
plot(FvCB_J_x, FvCB_J_y, 'Marker','square', ...   % solid line with square markers
                     'MarkerFaceColor', 'k', ... % k = black markers and line
                     'MarkerEdgeColor', 'k', ...
                     'MarkerSize', 8, ...
                     'Color','k',...
                     'DisplayName', 'FvCB');
hold on

% Plot Cc (x) against Gross Assimilation (y) for optimised alpha values (red dashed lines with red triangle markers)
plot(ePhoto_opt_x_J, ePhoto_opt_y_J, ...
            '--^', ...
            'Color', 'r', ...
            'MarkerFaceColor', 'r', ... % r = red markers and line
            'MarkerEdgeColor', 'r', ...
            'MarkerSize', 8, ...
            'Color','r',...
            'DisplayName', 'α = 1.12');

hold on

% Plot background dashed lines for comparison
plot(ePhoto_J_x_06, ePhoto_J_y_06, ...
            '--', ...
            'Color', '#E97132', ...
            'DisplayName', 'α = 0.6');
hold on
plot(ePhoto_J_x_08, ePhoto_J_y_08, ...
            '--', ...
            'Color', '#CC00FF', ...
            'DisplayName', 'α = 0.8');
hold on
plot(ePhoto_J_x_10, ePhoto_J_y_10, ...
            '--', ...
            'Color', '#00B0F0', ...
            'DisplayName', 'α = 1.0');
hold on
plot(ePhoto_J_x_12, ePhoto_J_y_12, ...
            '--', ...
            'Color', '#00B050', ...
            'DisplayName', 'α = 1.2');
hold on
plot(ePhoto_J_x_14, ePhoto_J_y_14, ...
            '--', ...
            'Color', '#0070C0', ...
            'DisplayName', 'α = 1.4');

% Set x axis limits
xlim([200 700]);
ylim([0 50]);

% Add labels, title, and legend
xlabel('C_c (μmol mol^{-1})','FontSize',12);
ylabel('Gross Assimilation (μmol m^{-2} s^{-1})','FontSize',12');
title('Scaling Factors for RuBP-Regeneration (α_{Enzymes}) ','FontSize', 12');

% Add legend, including FvCB
legend('show', 'Location', 'east');

%% Create figure for Rubisco carboxylation scaling factors
figure;
% Specify Cc and Gross Assimilation data for FvCB model
FvCB_V_x = Farq_Matrix_V(1:5,5); % Cc - only first five rows needed as they repeat
FvCB_V_y = Farq_Matrix_V(1:5,4); % Gross A - only first five rows needed as they repeat

% Specify Cc and Gross Assimilation data for optimal scaling factor
% Round off scaling factor values by 2 d.p. to prevent mismatches associated with floating-point numbers
ePhoto_Matrix_V(:, 1) = round(ePhoto_Matrix_V(:, 1), 2);

% Find rows for scaling_factors in ePhoto by finding and assigning indices of ePhotoMatrix from Vcmax_result
opt_indices_V = find(ePhoto_Matrix_V(:,1)==a_Rubisco);
idx06_V = find(ePhoto_Matrix_V(:,1)==0.6);
idx08_V = find(ePhoto_Matrix_V(:,1)==0.8);
idx10_V = find(ePhoto_Matrix_V(:,1)==1.0);
idx12_V = find(ePhoto_Matrix_V(:,1)==1.2);
idx14_V = find(ePhoto_Matrix_V(:,1)==1.4);

ePhoto_opt_x_V = ePhoto_Matrix_V(opt_indices_V,5); % Cc
ePhoto_opt_y_V = ePhoto_Matrix_V(opt_indices_V,4); % Gross A

ePhoto_V_x_06 = ePhoto_Matrix_V(idx06_V,5);
ePhoto_V_y_06 = ePhoto_Matrix_V(idx06_V,4);

ePhoto_V_x_08 = ePhoto_Matrix_V(idx08_V,5);
ePhoto_V_y_08 = ePhoto_Matrix_V(idx08_V,4);

ePhoto_V_x_10 = ePhoto_Matrix_V(idx10_V,5);
ePhoto_V_y_10 = ePhoto_Matrix_V(idx10_V,4);

ePhoto_V_x_12 = ePhoto_Matrix_V(idx12_V,5);
ePhoto_V_y_12 = ePhoto_Matrix_V(idx12_V,4);

ePhoto_V_x_14 = ePhoto_Matrix_V(idx14_V,5);
ePhoto_V_y_14 = ePhoto_Matrix_V(idx14_V,4);

% Plot FvCB points with markers and dashed lines
plot(FvCB_V_x, FvCB_V_y, 'Marker','square', ...   % solid line with square markers
                     'MarkerFaceColor', 'k', ... % k = black markers and line
                     'MarkerEdgeColor', 'k', ...
                     'MarkerSize', 8, ...
                     'Color','k',...
                     'DisplayName', 'FvCB');
hold on

% Plot Cc (x) against Gross Assimilation (y) for optimised alpha values (red dashed lines with red triangle markers)
plot(ePhoto_opt_x_V, ePhoto_opt_y_V, ...
            '--^', ...
            'Color', 'r', ...
            'MarkerFaceColor', 'r', ... % r = red markers and line
            'MarkerEdgeColor', 'r', ...
            'MarkerSize', 8, ...
            'Color','r',...
            'DisplayName', 'α = 1.32');

hold on

% Plot background dashed lines for comparison
plot(ePhoto_V_x_06, ePhoto_V_y_06, ...
            '--', ...
            'Color', '#E97132', ...
            'DisplayName', 'α = 0.6');
hold on
plot(ePhoto_V_x_08, ePhoto_V_y_08, ...
            '--', ...
            'Color', '#CC00FF', ...
            'DisplayName', 'α = 0.8');
hold on
plot(ePhoto_V_x_10, ePhoto_V_y_10, ...
            '--', ...
            'Color', '#00B0F0', ...
            'DisplayName', 'α = 1.0');
hold on
plot(ePhoto_V_x_12, ePhoto_V_y_12, ...
            '--', ...
            'Color', '#00B050', ...
            'DisplayName', 'α = 1.2');
hold on
plot(ePhoto_V_x_14, ePhoto_V_y_14, ...
            '--', ...
            'Color', '#0070C0', ...
            'DisplayName', 'α = 1.4');

% Set x axis limits
xlim([60 180]);
ylim([0 50]);

% Add labels, title, and legend
xlabel('C_c (μmol mol^{-1})','FontSize',12);
ylabel('Gross Assimilation (μmol m^{-2} s^{-1})','FontSize',12');
title('Scaling Factors for Rubisco Carboxylation (α_{Rubisco}) ','FontSize', 12');

% Add legend, including FvCB
legend('show', 'Location', 'northeast');