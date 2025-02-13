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
%% Create the violin plot for the different enzyme combinations
% figure;
% violin_plot = violinplot(CO2_data,'ShowData',false,'ShowMean',true);

% Adjust colors with triplet codes
%violinColors = ["#D95319";"#4DBEEE";"#77AC30";"#7E2F8E";"#F37581";"#F9BC76";"#F7F779";"#5BF17C";"#63DDA0";"#65B8F8"];
violinColors = [0.9686 0.9686 0.4745;...yellow - Stress at 130
                0.9765 0.7373 0.4627;...light orange - Current at 130
                0.9529 0.4588 0.5059;...red - Future at 130
                0.3882 0.8667 0.6275;...green - Stress at 250
                0.3569 0.9451 0.4863;...green - Current at 250
                0.4667 0.6745 0.1882;...green - Future at 250
                0.3961 0.7216 0.9725;...blue - Stress at 360
                0.8000 0.6000 1.0000;...purple - Current at 360
                0.4941 0.1843 0.5569];%purple - Future at 360
% Alternative - match colours used for strategies (Stress - red; current - magenta; future - blue)
% violinColors = [0.9529 0.4588 0.5059;... red - Stress at 130
%                 0.85 0.5 0.85;... magenta - Current at 130
%                 0.3961 0.7216 0.9725;... blue - Future at 130
%                 0.9529 0.4588 0.5059;... red - Stress at 250
%                 0.85 0.5 0.85;... magenta - Current at 250
%                 0.3961 0.7216 0.9725;... blue - Future at 250
%                 0.9529 0.4588 0.5059;... red - Stress at 360
%                 0.85 0.5 0.85;... magenta - Current at 360
%                 0.3961 0.7216 0.9725];%  blue- Future at 360

figure;
%hold on;
% Store handles for violin plots
violinHandles = cell(1, 9); 

% Loop through each group and create violin plots
for k = 1:9
    violinHandles{k} = Violin(CO2_data(:, k), k, 'ViolinColor', violinColors(k, :), 'ShowMean', true, 'ShowData', false);
end
%hold off;

% Add simple labels for the enzymes
set(gca, 'XTick', 1:9, 'XTickLabel', {"Stress",... %Rubisco SBPase 130
                                      "Current",... %Rubisco Aldolase SBPase PRK 130
                                      "Future",... %Rubisco Aldolase SBPase PRK FBPase TK 130
                                      "Stress",... %Rubisco SBPase 250
                                      "Current",... %Rubisco Aldolase SBPase PRK 250
                                      "Future",... %Rubisco Aldolase SBPase PRK FBPase TK 250
                                      "Stress",... %Rubisco SBPase 360
                                      "Current",... %Rubisco Aldolase SBPase PRK 360
                                      "Future"},'FontSize', 20); %Rubisco Aldolase SBPase PRK FBPase TK 360

groupLabels = {'C_{c} = 130 μmol mol^{−1}', 'C_{c} = 250 μmol mol^{−1}', 'C_{c} = 360 μmol mol^{−1}'};

% Manually position group labels above the violins
for i = 1:length(groupLabels)
    % Calculate the x-position as the mean of the x-coordinates of the violins in each group
    xPos = mean([violinHandles{(i-1)*3+1}.ScatterPlot.XData, ...
                 violinHandles{(i-1)*3+2}.ScatterPlot.XData, ...
                 violinHandles{(i-1)*3+3}.ScatterPlot.XData]);

    % Calculate the y-position slightly above the maximum y-value of the violins in each group
    % yPos = max([violinHandles{(i-1)*3+1}.ScatterPlot.YData, ...
    %             violinHandles{(i-1)*3+2}.ScatterPlot.YData, ...
    %             violinHandles{(i-1)*3+3}.ScatterPlot.YData]);
    yPos = 12;

    % Place the group label above the violins
    text(xPos, yPos, groupLabels{i}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 18);
end

% Add dashed horizontal line at y = 0
yline(0, '--k', 'LineWidth', 1.5);

% Rotate labels
xtickangle(45);

% Customize the labels and title
xlabel('Strategy','FontSize', 20);
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})','FontSize', 20); % absolute changes

% Specify limits for y axis
ylim([-5, 8]); % absolute changes

% Find mean, min and max of percentages
mean_percents = mean(CO2_data_percent);
min_max_percents = vertcat(min(CO2_data_percent),max(CO2_data_percent));

% 
% for k = 1:11
%     h = findobj(gca, 'Type', 'Patch');
%     set(h(k), 'FaceColor', violinColors(k, :));
% end

% Save and print the plot
set(gcf, 'PaperOrientation', 'landscape');
%print(gcf,fullfile('Outputs/rice_params/graphs',"violin_enzymes_new3.pdf"),'-dpdf','-bestfit');

% Number of percentages greater than or equal to 0
% pos_percents = sum(CO2_data_percent >= 0);
% neg_percents = 2000-pos_percents;

%Plot stacked bar to view proportion of positive and negative percentage changes
%prop_percents = vertcat(pos_percents,neg_percents)';
%figure;
%bar(categorical({'RS 130','RS 250','RS 360','RSAP 130','RSAP 250','RSAP 360','RSAPFT 130','RSAPFT 250','RSAPFT 360'}),prop_percents,'stacked');
%legend('>0', '<0','Location','eastoutside');
%title('Proportion of Positive/Negative Percentage Changes in Gross Assimilation')
%xlabel('Enzymes');
%ylabel('No. of Fold Changes');
%set(gcf, 'PaperOrientation', 'landscape');
%print(gcf,fullfile('Outputs/rice_params/graphs',"prop_enzymes_percent.pdf"),'-dpdf','-bestfit');