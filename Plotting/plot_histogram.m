% Produce a histogram for enzyme sensitivity analysis 
%% Load data
load output_enzyme_adjustment_test_250_2000_3.mat; % Load results from 2000 node test simulation
%% Create outputs for increases in gross CO2 assimilation rate
[IncreaseInGrossAssimilationRate_RSAP,... %four
    IncreaseInGrossAssimilationRate_RA,IncreaseInGrossAssimilationRate_RS,IncreaseInGrossAssimilationRate_RP,...
        IncreaseInGrossAssimilationRate_AS,IncreaseInGrossAssimilationRate_AP,IncreaseInGrossAssimilationRate_SP,...%two
        IncreaseInGrossAssimilationRate_R,IncreaseInGrossAssimilationRate_A,IncreaseInGrossAssimilationRate_S,IncreaseInGrossAssimilationRate_P] = deal(zeros(2000,1));
%% Convert values for plotting
%Expressing the y-axis as the increases/decreases in gross CO2 assimilation rate resulting from FC = 1  (first row in each vector of assimilation rates)
for n=1:length(GrossAssimilationRate_RSAP)
IncreaseInGrossAssimilationRate_RSAP(n) = GrossAssimilationRate_RSAP(n)-GrossAssimilationRate_RSAP(1);
IncreaseInGrossAssimilationRate_RA(n) = GrossAssimilationRate_RA(n)-GrossAssimilationRate_RA(1);
IncreaseInGrossAssimilationRate_RS(n) = GrossAssimilationRate_RS(n)-GrossAssimilationRate_RS(1);
IncreaseInGrossAssimilationRate_RP(n) = GrossAssimilationRate_RP(n)-GrossAssimilationRate_RP(1);
IncreaseInGrossAssimilationRate_AS(n) = GrossAssimilationRate_AS(n)-GrossAssimilationRate_AS(1);
IncreaseInGrossAssimilationRate_AP(n) = GrossAssimilationRate_AP(n)-GrossAssimilationRate_AP(1);
IncreaseInGrossAssimilationRate_SP(n) = GrossAssimilationRate_SP(n)-GrossAssimilationRate_SP(1);
IncreaseInGrossAssimilationRate_R(n) = GrossAssimilationRate_R(n)-GrossAssimilationRate_R(1);
IncreaseInGrossAssimilationRate_A(n) = GrossAssimilationRate_A(n)-GrossAssimilationRate_A(1);
IncreaseInGrossAssimilationRate_S(n) = GrossAssimilationRate_S(n)-GrossAssimilationRate_S(1);
IncreaseInGrossAssimilationRate_P(n) = GrossAssimilationRate_P(n)-GrossAssimilationRate_P(1);
end

%Expressing the x-axis as the multiplicative change in protein content - difference in fold change from non optimal to optimal
change_aldolase_FC = aldolase_FC-1;
change_rubisco_FC = rubisco_FC-1;
change_sbpase_FC = sbpase_FC-1;
change_prk_FC = prk_FC-1;

% Get the min and max values for changes in enzyme FC
plotting_range_FC_R = horzcat(min(change_rubisco_FC),max(change_rubisco_FC));
plotting_range_FC_A = horzcat(min(change_aldolase_FC),max(change_aldolase_FC));
plotting_range_FC_S = horzcat(min(change_sbpase_FC),max(change_sbpase_FC));
plotting_range_FC_P = horzcat(min(change_prk_FC),max(change_prk_FC));

% Get the min and max values for changes in assimilation rates 
plotting_range_R = horzcat(min(IncreaseInGrossAssimilationRate_R),max(IncreaseInGrossAssimilationRate_R));
plotting_range_A = horzcat(min(IncreaseInGrossAssimilationRate_A),max(IncreaseInGrossAssimilationRate_A));
plotting_range_S = horzcat(min(IncreaseInGrossAssimilationRate_S),max(IncreaseInGrossAssimilationRate_S));
plotting_range_P = horzcat(min(IncreaseInGrossAssimilationRate_P),max(IncreaseInGrossAssimilationRate_P));
plotting_range_RA = horzcat(min(IncreaseInGrossAssimilationRate_RA),max(IncreaseInGrossAssimilationRate_RA));
plotting_range_RS = horzcat(min(IncreaseInGrossAssimilationRate_RS),max(IncreaseInGrossAssimilationRate_RS));
plotting_range_RP = horzcat(min(IncreaseInGrossAssimilationRate_RP),max(IncreaseInGrossAssimilationRate_RP));
plotting_range_AS = horzcat(min(IncreaseInGrossAssimilationRate_AS),max(IncreaseInGrossAssimilationRate_AS));
plotting_range_AP = horzcat(min(IncreaseInGrossAssimilationRate_AP),max(IncreaseInGrossAssimilationRate_AP));
plotting_range_SP = horzcat(min(IncreaseInGrossAssimilationRate_SP),max(IncreaseInGrossAssimilationRate_SP));
plotting_range_RSAP = horzcat(min(IncreaseInGrossAssimilationRate_RSAP),max(IncreaseInGrossAssimilationRate_RSAP));

%% Prepare data for histogram
% Create a matrix with all the increases/decreases in assimilation rates
CO2_data = horzcat(IncreaseInGrossAssimilationRate_R,...
                 IncreaseInGrossAssimilationRate_A,...
                 IncreaseInGrossAssimilationRate_S,...
                 IncreaseInGrossAssimilationRate_P,...
                 IncreaseInGrossAssimilationRate_RA,...
                 IncreaseInGrossAssimilationRate_RS,...
                 IncreaseInGrossAssimilationRate_RP,...
                 IncreaseInGrossAssimilationRate_AS,...
                 IncreaseInGrossAssimilationRate_AP,...
                 IncreaseInGrossAssimilationRate_SP,...
                 IncreaseInGrossAssimilationRate_RSAP);
%% Create the histogram for the different enzyme combinations
enzyme_names = string({"Rubisco","Aldolase","SBPase","PRK",...
                       "Rubisco Aldolase","Rubisco SBPase",...
                       "Rubisco PRK","Aldolase SBPase",...
                       "Aldolase PRK","SBPase PRK","Rubisco Aldolase SBPase PRK"});
enzyme_colours = string({"#D95319","#4DBEEE","#77AC30","#7E2F8E", ... single enzyme colours
                        "#F37581","#F9BC76","#F7F779","#5BF17C","#63DDA0","#65B8F8"}); % paired enzyme colours

titleName = strcat ("Distribution of 2000 FC for",enzyme_names);

% Plot single enzyme distributions
figure
tiledlayout(2,2);
for i = 1:4 % First set of indices for assimilation rate changes
    nexttile
    histogram(CO2_data(:,i),'BinEdges',-4:0.5:10,'FaceColor',enzyme_colours(i))
    title(strcat ("Distribution of 2000 FC for ",enzyme_names(i)))
    xlabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})')
    ylabel('Frequency')
end
% Save the plot
set(gcf, 'PaperOrientation', 'landscape');
print(gcf,fullfile('Outputs/rice_params/graphs',"histogram_single_enzymes.pdf"),'-dpdf','-fillpage');

% Plot paired enzyme distributions
figure
tiledlayout(3,3);
for j = 5:10 % Second set of indices for assimilation rate changes
    nexttile
    histogram(CO2_data(:,j),'BinEdges',-4:0.5:10,'FaceColor',enzyme_colours(j))
    title(strcat ("Distribution of 2000 FC for ",enzyme_names(j)))
    xlabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})')
    ylabel('Frequency')
end
% Save the plot
set(gcf, 'PaperOrientation', 'landscape');
print(gcf,fullfile('Outputs/rice_params/graphs',"histogram_paired_enzymes.pdf"),'-dpdf','-fillpage');

% Plot four enzyme distribution
figure
histogram(CO2_data(:,11),'BinEdges',-4:0.5:10,'FaceColor',"#CC99FF")
title(strcat ("Distribution of 2000 FC for ",enzyme_names(11)))
xlabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})')
ylabel('Frequency')
% Save the plot
set(gcf, 'PaperOrientation', 'landscape');
print(gcf,fullfile('Outputs/rice_params/graphs',"histogram_quad_enzymes.pdf"),'-dpdf','-fillpage');
