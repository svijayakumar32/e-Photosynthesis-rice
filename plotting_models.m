% Plot results of plot_models
load('Farq_ePhoto_models.mat')

% Create figure
figure;
hold on;

% Specify Cc and Gross Assimilation data for FvCB model
FvCB_x = Farq_Matrix(:,5); %Cc
FvCB_y = Farq_Matrix(:,6); %Gross A

% Plot FvCB points with markers and dashed lines
plot(FvCB_x, FvCB_y, 'Marker','square', ...   % '--o' specifies dashed line with circle markers
                     'MarkerFaceColor', 'k', ... % k=black markers and line
                     'MarkerEdgeColor', 'k', ...
                     'MarkerSize', 8, ...
                     'Color','k',...
                     'DisplayName', 'FvCB');

% Specify Cc and Gross Assimilation data for non-tuned ePhotosynthesis
% model (both alphas = 1)
ePhoto_x = ePhoto_Matrix_1(:,5);
ePhoto_y = ePhoto_Matrix_1(:,6);

ePhoto_opt_x = ePhoto_Matrix_2(:,5);
ePhoto_opt_y = ePhoto_Matrix_2(:,6);

% Plot Cc (x) against Gross Assimilation (y) for alpha values=1 (black dashed lines)
plot(ePhoto_x, ePhoto_y, ...
            '--o', ...
            'Color', 'k', ...
            'MarkerFaceColor', 'k', ... % k=black markers and line
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 8, ...
            'Color','k',...
            'DisplayName', 'ePhotosynthesis (α_Enzymes = 1, α_Rubisco = 1)');

% Plot Cc (x) against Gross Assimilation (y) for optimised alpha values (red dashed lines)
plot(ePhoto_opt_x, ePhoto_opt_y, ...
            '--o', ...
            'Color', 'r', ...
            'MarkerFaceColor', 'r', ... % k=black markers and line
            'MarkerEdgeColor', 'r', ...
            'MarkerSize', 8, ...
            'Color','r',...
            'DisplayName', 'ePhotosynthesis (α_Enzymes = 0.98, α_Rubisco = 0.96)');
hold off;
% 
% % Set the x-axis limits from 0 to 700
xlim([0 1000]);
set(gca, 'FontSize', 20); 
% 
% % Add labels, title, and legend
xlabel('C_c (μmol mol^{-1})','FontSize',20);
ylabel('Gross Assimilation (μmol m^{-2} s^{-1})','FontSize',20');
title('Scaling Factors for Rubisco and RuBP-Regeneration Limitation','FontSize',12');

% Add legend, including FvCB
legend('show', 'Location', 'east');
