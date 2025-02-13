% Produce violin plots for enzyme sensitivity analysis using the violinplot
% function from cobratoolbox (https://github.com/opencobra/cobratoolbox/blob/master/external/base/utilities/violinPlots/violinplot.m)
%% Load data
load output_enzyme_adjustment_test_2000_5_new.mat; % Load results from 2000 node test simulation
%% Create outputs for increases in gross CO2 assimilation rate
[IncreaseInGrossAssimilationRate_RS_130,...%two
 IncreaseInGrossAssimilationRate_RS_250,...%two
 IncreaseInGrossAssimilationRate_RS_360,...%two
    IncreaseInGrossAssimilationRate_RSAP_130,... %four
    IncreaseInGrossAssimilationRate_RSAP_250,... %four
    IncreaseInGrossAssimilationRate_RSAP_360,... %four
        IncreaseInGrossAssimilationRate_RSAPFT_130,...
        IncreaseInGrossAssimilationRate_RSAPFT_250,...
        IncreaseInGrossAssimilationRate_RSAPFT_360] = deal(zeros(2000,1));
%% To express this as % change
[PercentChangeGrossAssimilationRate_RS_130,...%two
 PercentChangeGrossAssimilationRate_RS_250,...%two
 PercentChangeGrossAssimilationRate_RS_360,...%two
    PercentChangeGrossAssimilationRate_RSAP_130,... %four
    PercentChangeGrossAssimilationRate_RSAP_250,... %four
    PercentChangeGrossAssimilationRate_RSAP_360,... %four
        PercentChangeGrossAssimilationRate_RSAPFT_130,...
        PercentChangeGrossAssimilationRate_RSAPFT_250,...
        PercentChangeGrossAssimilationRate_RSAPFT_360] = deal(zeros(2000,1));
%% Convert values for plotting
%Expressing the y-axis as the increases/decreases in gross CO2 assimilation rate resulting from FC = 1  (first row in each vector of assimilation rates)
for n=1:length(GrossAssimilationRate_RS_130)
    IncreaseInGrossAssimilationRate_RS_130(n) = GrossAssimilationRate_RS_130(n)-GrossAssimilationRate_RS_130(1);
    IncreaseInGrossAssimilationRate_RS_250(n) = GrossAssimilationRate_RS_250(n)-GrossAssimilationRate_RS_250(1);
    IncreaseInGrossAssimilationRate_RS_360(n) = GrossAssimilationRate_RS_360(n)-GrossAssimilationRate_RS_360(1);
    IncreaseInGrossAssimilationRate_RSAP_130(n) = GrossAssimilationRate_RSAP_130(n)-GrossAssimilationRate_RSAP_130(1);
    IncreaseInGrossAssimilationRate_RSAP_250(n) = GrossAssimilationRate_RSAP_250(n)-GrossAssimilationRate_RSAP_250(1);
    IncreaseInGrossAssimilationRate_RSAP_360(n) = GrossAssimilationRate_RSAP_360(n)-GrossAssimilationRate_RSAP_360(1);
    IncreaseInGrossAssimilationRate_RSAPFT_130(n) = GrossAssimilationRate_RSAPFT_130(n)-GrossAssimilationRate_RSAPFT_130(1);
    IncreaseInGrossAssimilationRate_RSAPFT_250(n) = GrossAssimilationRate_RSAPFT_250(n)-GrossAssimilationRate_RSAPFT_250(1);
    IncreaseInGrossAssimilationRate_RSAPFT_360(n) = GrossAssimilationRate_RSAPFT_360(n)-GrossAssimilationRate_RSAPFT_360(1);
end
%OR as percent change
for n=1:length(GrossAssimilationRate_RS_130)
    PercentChangeGrossAssimilationRate_RS_130(n) = (GrossAssimilationRate_RS_130(n)-GrossAssimilationRate_RS_130(1))/GrossAssimilationRate_RS_130(1)*100;
    PercentChangeGrossAssimilationRate_RS_250(n) = (GrossAssimilationRate_RS_250(n)-GrossAssimilationRate_RS_250(1))/GrossAssimilationRate_RS_250(1)*100;
    PercentChangeGrossAssimilationRate_RS_360(n) = (GrossAssimilationRate_RS_360(n)-GrossAssimilationRate_RS_360(1))/GrossAssimilationRate_RS_360(1)*100;
    PercentChangeGrossAssimilationRate_RSAP_130(n) = (GrossAssimilationRate_RSAP_130(n)-GrossAssimilationRate_RSAP_130(1))/GrossAssimilationRate_RSAP_130(1)*100;
    PercentChangeGrossAssimilationRate_RSAP_250(n) = (GrossAssimilationRate_RSAP_250(n)-GrossAssimilationRate_RSAP_250(1))/GrossAssimilationRate_RSAP_250(1)*100;
    PercentChangeGrossAssimilationRate_RSAP_360(n) = (GrossAssimilationRate_RSAP_360(n)-GrossAssimilationRate_RSAP_360(1))/GrossAssimilationRate_RSAP_360(1)*100;
    PercentChangeGrossAssimilationRate_RSAPFT_130(n) = (GrossAssimilationRate_RSAPFT_130(n)-GrossAssimilationRate_RSAPFT_130(1))/GrossAssimilationRate_RSAPFT_130(1)*100;
    PercentChangeGrossAssimilationRate_RSAPFT_250(n) = (GrossAssimilationRate_RSAPFT_250(n)-GrossAssimilationRate_RSAPFT_250(1))/GrossAssimilationRate_RSAPFT_250(1)*100;
    PercentChangeGrossAssimilationRate_RSAPFT_360(n) = (GrossAssimilationRate_RSAPFT_360(n)-GrossAssimilationRate_RSAPFT_360(1))/GrossAssimilationRate_RSAPFT_360(1)*100;
end

%Expressing the x-axis as the multiplicative change in protein content - difference in fold change from initial 
change_RS_rubisco_FC = RS_rubisco_FC-1;
change_RS_sbpase_FC = RS_sbpase_FC-1;
change_RSAP_rubisco_FC = RSAP_rubisco_FC-1;
change_RSAP_sbpase_FC = RSAP_sbpase_FC-1;
change_RSAP_aldolase_FC = RSAP_aldolase_FC-1;
change_RSAP_prk_FC = RSAP_prk_FC-1;
change_RSAPFT_rubisco_FC = RSAPFT_rubisco_FC-1;
change_RSAPFT_sbpase_FC = RSAPFT_sbpase_FC-1;
change_RSAPFT_aldolase_FC = RSAPFT_aldolase_FC-1;
change_RSAPFT_prk_FC = RSAPFT_prk_FC-1;
change_RSAPFT_fbpase_FC = RSAPFT_fbpase_FC-1;
change_RSAPFT_tk_FC = RSAPFT_tk_FC-1;

% Get the min and max values for changes in enzyme FC
plotting_range_FC_RS_rubisco = horzcat(min(change_RS_rubisco_FC),max(change_RS_rubisco_FC));
plotting_range_FC_RS_sbpase = horzcat(min(change_RS_sbpase_FC),max(change_RS_rubisco_FC));
plotting_range_FC_RSAP_rubisco = horzcat(min(change_RSAP_rubisco_FC),max(change_RSAP_rubisco_FC));
plotting_range_FC_RSAP_sbpase = horzcat(min(change_RSAP_sbpase_FC),max(change_RSAP_sbpase_FC));
plotting_range_FC_RSAP_aldolase = horzcat(min(change_RSAP_aldolase_FC),max(change_RSAP_aldolase_FC));
plotting_range_FC_RSAP_prk = horzcat(min(change_RSAP_prk_FC),max(change_RSAP_prk_FC));

plotting_range_FC_RSAPFT_rubisco = horzcat(min(change_RSAPFT_rubisco_FC),max(change_RSAPFT_rubisco_FC));
plotting_range_FC_RSAPFT_sbpase = horzcat(min(change_RSAPFT_sbpase_FC),max(change_RSAPFT_sbpase_FC));
plotting_range_FC_RSAPFT_aldolase = horzcat(min(change_RSAPFT_aldolase_FC),max(change_RSAPFT_aldolase_FC));
plotting_range_FC_RSAPFT_prk = horzcat(min(change_RSAPFT_prk_FC),max(change_RSAPFT_prk_FC));
plotting_range_FC_RSAPFT_fbpase = horzcat(min(change_RSAPFT_fbpase_FC),max(change_RSAPFT_fbpase_FC));
plotting_range_FC_RSAPFT_tk = horzcat(min(change_RSAPFT_tk_FC),max(change_RSAPFT_tk_FC));

% Get the min and max values for changes in assimilation rates 
plotting_range_RS_130 = horzcat(min(IncreaseInGrossAssimilationRate_RS_130),max(IncreaseInGrossAssimilationRate_RS_130));
plotting_range_RS_250 = horzcat(min(IncreaseInGrossAssimilationRate_RS_250),max(IncreaseInGrossAssimilationRate_RS_250));
plotting_range_RS_360 = horzcat(min(IncreaseInGrossAssimilationRate_RS_360),max(IncreaseInGrossAssimilationRate_RS_360));
plotting_range_RSAP_130 = horzcat(min(IncreaseInGrossAssimilationRate_RSAP_130),max(IncreaseInGrossAssimilationRate_RSAP_130));
plotting_range_RSAP_250 = horzcat(min(IncreaseInGrossAssimilationRate_RSAP_250),max(IncreaseInGrossAssimilationRate_RSAP_250));
plotting_range_RSAP_360 = horzcat(min(IncreaseInGrossAssimilationRate_RSAP_360),max(IncreaseInGrossAssimilationRate_RSAP_360));
plotting_range_RSAPFT_130 = horzcat(min(IncreaseInGrossAssimilationRate_RSAPFT_130),max(IncreaseInGrossAssimilationRate_RSAPFT_130));
plotting_range_RSAPFT_250 = horzcat(min(IncreaseInGrossAssimilationRate_RSAPFT_250),max(IncreaseInGrossAssimilationRate_RSAPFT_250));
plotting_range_RSAPFT_360 = horzcat(min(IncreaseInGrossAssimilationRate_RSAPFT_360),max(IncreaseInGrossAssimilationRate_RSAPFT_360));

%% Prepare data for violin plot
% Create a table with all the increases/decreases in assimilation rates
CO2_data = horzcat(IncreaseInGrossAssimilationRate_RS_130,...
                   IncreaseInGrossAssimilationRate_RS_250,...                   
                   IncreaseInGrossAssimilationRate_RS_360,...
                   IncreaseInGrossAssimilationRate_RSAP_130,...
                   IncreaseInGrossAssimilationRate_RSAP_250,...
                   IncreaseInGrossAssimilationRate_RSAP_360,...
                   IncreaseInGrossAssimilationRate_RSAPFT_130,...
                   IncreaseInGrossAssimilationRate_RSAPFT_250,...
                   IncreaseInGrossAssimilationRate_RSAPFT_360);
% OR as assimilation rates as percent changes
CO2_data_percent = horzcat(PercentChangeGrossAssimilationRate_RS_130,...
                   PercentChangeGrossAssimilationRate_RS_250,...
                   PercentChangeGrossAssimilationRate_RS_360,...
                   PercentChangeGrossAssimilationRate_RSAP_130,...
                   PercentChangeGrossAssimilationRate_RSAP_250,...
                   PercentChangeGrossAssimilationRate_RSAP_360,...
                   PercentChangeGrossAssimilationRate_RSAPFT_130,...
                   PercentChangeGrossAssimilationRate_RSAPFT_250,...
                   PercentChangeGrossAssimilationRate_RSAPFT_360);
% Find mean, min and max of percentages
mean_percents = mean(CO2_data_percent);
min_max_percents = vertcat(min(CO2_data_percent),max(CO2_data_percent));
%% Create the violin plot for the different enzyme combinations
% Set colors with triplet codes
violinColors_Stress = [0.9529 0.4588 0.5059; %Red
                0.9529 0.4588 0.5059; 
                0.9529 0.4588 0.5059]; 

violinColors_Current = [0.85 0.5 0.85; %Magenta
                 0.85 0.5 0.85; 
                 0.85 0.5 0.85]; 

violinColors_Future = [0.3961 0.7216 0.9725; % Blue
                0.3961 0.7216 0.9725;
                0.3961 0.7216 0.9725]; 

groupLabels = {'130', '250', '360'};

%% STRESS
figure;
for i = 1:3
    Violin(CO2_data(:, (i-1)*3+1), i, 'ViolinColor', violinColors_Stress(i, :), 'ShowMean', true, 'ShowData', false);
end
set(gca, 'XTick', 1:3, 'XTickLabel', groupLabels, 'FontSize', 20);
% Customize the labels and title
title('Stress', 'FontSize', 22);
xlabel('C_{c} (μmol mol^{−1})', 'FontSize', 20);
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})', 'FontSize', 20);
% Add dashed horizontal line at y = 0
yline(0, '--k', 'LineWidth', 1.5);
% Specify limits for y axis
ylim([-3, 6]);
% Rotate labels
xtickangle(45);
% Save and print the plot
set(gcf, 'PaperOrientation', 'landscape');
print(gcf,fullfile('Outputs/rice_params/graphs',"violin_enzymes_Stress.pdf"),'-dpdf','-bestfit');figure;

%% CURRENT
for i = 1:3
    Violin(CO2_data(:, (i-1)*3+2), i, 'ViolinColor', violinColors_Current(i, :), 'ShowMean', true, 'ShowData', false);
end
set(gca, 'XTick', 1:3, 'XTickLabel', groupLabels, 'FontSize', 20);
% Customize the labels and title
title('Current', 'FontSize', 22);
xlabel('C_{c} (μmol mol^{−1})', 'FontSize', 20);
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})', 'FontSize', 20);
% Add dashed horizontal line at y = 0
yline(0, '--k', 'LineWidth', 1.5);
% Specify limits for y axis
ylim([-3, 6]);
% Rotate labels
xtickangle(45);
% Save and print the plot
set(gcf, 'PaperOrientation', 'landscape');
print(gcf,fullfile('Outputs/rice_params/graphs',"violin_enzymes_Current.pdf"),'-dpdf','-bestfit');figure;

%% FUTURE
for i = 1:3
    Violin(CO2_data(:, (i-1)*3+3), i, 'ViolinColor', violinColors_Future(i, :), 'ShowMean', true, 'ShowData', false);
end
set(gca, 'XTick', 1:3, 'XTickLabel', groupLabels, 'FontSize', 20);
% Customize the labels and title
title('Future', 'FontSize', 22);
xlabel('C_{c} (μmol mol^{−1})', 'FontSize', 20);
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})', 'FontSize', 20);
% Add dashed horizontal line at y = 0
yline(0, '--k', 'LineWidth', 1.5);
% Specify limits for y axis
ylim([-3, 6]);
% Rotate labels
xtickangle(45);
% Save and print the plot
set(gcf, 'PaperOrientation', 'landscape');
print(gcf,fullfile('Outputs/rice_params/graphs',"violin_enzymes_Future.pdf"),'-dpdf','-bestfit');

