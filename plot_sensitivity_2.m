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
plotting_range_R = horzcat(min(change_rubisco_FC),max(change_rubisco_FC));
plotting_range_A = horzcat(min(change_aldolase_FC),max(change_aldolase_FC));
plotting_range_S = horzcat(min(change_sbpase_FC),max(change_sbpase_FC));
plotting_range_P = horzcat(min(change_prk_FC),max(change_prk_FC));

% Get the min and max values for changes in assimilation rates
plotting_range_RA = horzcat(min(IncreaseInGrossAssimilationRate_RA),max(IncreaseInGrossAssimilationRate_RA));
plotting_range_RS = horzcat(min(IncreaseInGrossAssimilationRate_RS),max(IncreaseInGrossAssimilationRate_RS));
plotting_range_RP = horzcat(min(IncreaseInGrossAssimilationRate_RP),max(IncreaseInGrossAssimilationRate_RP));
plotting_range_AP = horzcat(min(IncreaseInGrossAssimilationRate_AP),max(IncreaseInGrossAssimilationRate_AP));
plotting_range_AS = horzcat(min(IncreaseInGrossAssimilationRate_AS),max(IncreaseInGrossAssimilationRate_AS));
plotting_range_SP = horzcat(min(IncreaseInGrossAssimilationRate_SP),max(IncreaseInGrossAssimilationRate_SP));
plotting_range_RSAP = horzcat(min(IncreaseInGrossAssimilationRate_RSAP),max(IncreaseInGrossAssimilationRate_RSAP));

% Options to plot two-enzyme changes and the resulting assimilation rates are 3D scatter or stem graphs, e.g. 

% Plotting a small set of dummy values
% X = [1,2,5,7,4]; %Rubisco
% Y = [3,2,5,3,4]; %Aldolase
% Z = [23,20,26,29,24]; %GrossAssimilation_RA
% figure
% stem3(X,Y,Z)

%% Scatter graph with a full set of 2000 nodes for Rubisco only
%figure
o = tiledlayout(2,2,'Padding','compact');
title(o,'Single enzyme increases vs. changes in assimilation rate')
nexttile
scatter(change_rubisco_FC,IncreaseInGrossAssimilationRate_R,"filled");
% Plot optimal point on top for comparison
hold on 
scatter(change_rubisco_FC(2), IncreaseInGrossAssimilationRate_R(2), 150,"pentagram",'filled','MarkerFaceColor',"#D95319");% plot optimal Rubisco on top in red
xticks([0 0.2 0.4 0.6 0.8 1]);
yticks([-3 -2 -1 0 1 2 3]);
xlabel('Additional Rubisco');
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
axis padded;

% Scatter graph with a full set of 2000 nodes for Aldolase only
nexttile
scatter(change_aldolase_FC,IncreaseInGrossAssimilationRate_A,"filled");
% Plot optimal point on top for comparison
hold on 
scatter(change_aldolase_FC(2), IncreaseInGrossAssimilationRate_A(2), 150,"pentagram",'filled','MarkerFaceColor',"#EDB120");% plot optimal aldolase on top in yellow
xticks([0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4]);
%yticks([-3 -2 -1 0 1 2 3]);
xlabel('Additional Aldolase');
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
axis padded;

% Scatter graph with a full set of 2000 nodes for SBPase only
nexttile
scatter(change_sbpase_FC,IncreaseInGrossAssimilationRate_S,"filled");
% Plot optimal point on top for comparison
hold on 
scatter(change_sbpase_FC(2), IncreaseInGrossAssimilationRate_S(2), 150,"pentagram",'filled','MarkerFaceColor',"#77AC30");% plot optimal SBPase on top in green
xticks([0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4 4.8 5.2 5.6]);
%yticks([-3 -2 -1 0 1 2 3]);
xlabel('Additional SBPase');
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
axis padded;

% Scatter graph with a full set of 2000 nodes for PRK only
nexttile
scatter(change_prk_FC,IncreaseInGrossAssimilationRate_P,"filled");
% Plot optimal point on top for comparison
hold on 
scatter(change_prk_FC(2), IncreaseInGrossAssimilationRate_P(2), 150,"pentagram",'filled','MarkerFaceColor',"#7E2F8E");% plot optimal PRK on top in purple
xticks([0 0.4 0.8 1.2 1.6 2.0 2.4 2.8]);
%yticks([-3 -2 -1 0 1 2 3]);
xlabel('Additional PRK');
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
axis padded;
set(gcf, 'PaperOrientation', 'landscape');
print(gcf,fullfile('Outputs/rice_params/graphs',"sensitivity_single_enzymes_new.pdf"),'-dpdf','-fillpage');
%% Scatter graph with a full set of 2000 nodes for Rubisco-Aldolase pair
figure
t = tiledlayout(3,3,'Padding','compact');
title(t,'Pairwise enzyme increases vs. changes in assimilation rate')
nexttile
% X is Rubisco, Y is Aldolase, Z is Assimilation
scatter3(change_rubisco_FC,change_aldolase_FC,IncreaseInGrossAssimilationRate_RA,"filled",'MarkerFaceAlpha','0.2','MarkerEdgeAlpha','0.2');
% Plot optimal point on top for comparison
hold on 
scatter3(change_rubisco_FC(2),change_aldolase_FC(2),IncreaseInGrossAssimilationRate_RA(2), 150,"pentagram",'filled','MarkerFaceColor',"#EDB120",'MarkerEdgeColor', 'black','MarkerFaceAlpha','1','MarkerEdgeAlpha','1');% plot optimal point on top in red
%scatter(change_aldolase_FC(2),IncreaseInGrossAssimilationRate_RA(2), 150,"pentagram",'filled','MarkerFaceColor',"#EDB120");% plot optimal Aldolase top in yellow
xticks([0 0.4 0.8 1.2]);
yticks([0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4]);
zticks([-3.5 -3.0 -2.5 -2.0 -1.5 -1.0 -0.5 0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5]);
xlabel('Additional Rubisco','Rotation',10);
ylabel('Additional Aldolase','Rotation',-20);
zlabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
axis padded;

%% Scatter graph with a full set of 2000 nodes for Rubisco-SBPase pair
nexttile 
% X is Rubisco, Y is SBPase, Z is Assimilation
scatter3(change_rubisco_FC,change_sbpase_FC,IncreaseInGrossAssimilationRate_RS,"filled",'MarkerFaceAlpha','0.2','MarkerEdgeAlpha','0.2'); % y = increase in GA or GA itself? ΔCO_{2} is a change
% Plot optimal point on top for comparison
hold on 
scatter3(change_rubisco_FC(2), change_sbpase_FC(2), IncreaseInGrossAssimilationRate_RS(2), 150,"pentagram",'filled','MarkerFaceColor',"#EDB120",'MarkerEdgeColor','black','MarkerFaceAlpha','1','MarkerEdgeAlpha','1');% plot optimal Rubisco on top in red 
%scatter(change_sbpase_FC(2), IncreaseInGrossAssimilationRate_RS(2), 150,"pentagram",'filled','MarkerFaceColor',"#77AC30");% plot optimal SBPase on top in green
xticks([0 0.2 0.4 0.6 0.8 1.0]);
yticks([0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4 4.8 5.2 5.6 6.0]);
zticks([-2.5 -2.0 -1.5 -1.0 -0.5 0 0.5 1.0 1.5 2.0 2.5 3.0 3.5]);
xlabel('Additional Rubisco','Rotation',10);
ylabel('Additional SBPase','Rotation',-20);
zlabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
axis padded;

%% Scatter graph with a full set of 2000 nodes for Rubisco-PRK pair
nexttile
% X is Rubisco, Y is PRK, Z is Assimilation
scatter3(change_rubisco_FC,change_prk_FC,IncreaseInGrossAssimilationRate_RP,"filled",'MarkerFaceAlpha','0.2','MarkerEdgeAlpha','0.2');
% Plot optimal point on top for comparison
hold on 
scatter3(change_rubisco_FC(2),change_prk_FC(2),IncreaseInGrossAssimilationRate_RP(2), 150,"pentagram",'filled','MarkerFaceColor',"#EDB120",'MarkerEdgeColor', 'black','MarkerFaceAlpha','1','MarkerEdgeAlpha','1');% plot optimal Rubisco on top in red
%scatter(change_prk_FC(2),IncreaseInGrossAssimilationRate_RP(2), 150,"pentagram",'filled','MarkerFaceColor',"#7E2F8E");% plot optimal PRK on top in purple
xticks([0 0.2 0.4 0.6 0.8 1.0]);
yticks([0 0.4 0.8 1.2 1.6 2.0 2.4 2.8]);
zticks([-3.5 -3.0 -2.5 -2.0 -1.5 -1.0 -0.5 0 0.5 1.0 1.5 2.0 2.5 3.0 3.5]);
xlabel('Additional Rubisco','Rotation',10);
ylabel('Additional PRK','Rotation',-20);
zlabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
axis padded;

%% Scatter graph with a full set of 2000 nodes for Aldolase-SBPase pair
nexttile
% X is Aldolase, Y is SBPase, Z is Assimilation
scatter3(change_aldolase_FC,change_sbpase_FC,IncreaseInGrossAssimilationRate_AS,"filled",'MarkerFaceAlpha','0.2','MarkerEdgeAlpha','0.2');
% Plot optimal point on top for comparison
hold on 
scatter3(change_aldolase_FC(2), change_sbpase_FC(2), IncreaseInGrossAssimilationRate_AS(2), 150,"pentagram",'filled','MarkerFaceColor',"#EDB120",'MarkerEdgeColor', 'black','MarkerFaceAlpha','1','MarkerEdgeAlpha','1');% plot optimal Aldolase on top in yellow
%scatter(change_sbpase_FC(2), IncreaseInGrossAssimilationRate_AS(2), 150,"pentagram",'filled','MarkerFaceColor',"#77AC30");% plot optimal SBPase on top in green
xticks([0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4]);
yticks([0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4 4.8 5.2 5.6 6.0]);
zticks([0 0.2 0.4 0.6 0.8 1.0 1.2]);
xlabel('Additional Aldolase','Rotation',10);
ylabel('Additional SBPase','Rotation',-20);
zlabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
axis padded;

%% Scatter graph with a full set of 2000 nodes for Aldolase-PRK pair
nexttile
% X is Aldolase, Y is PRK, Z is Assimilation
scatter3(change_aldolase_FC,change_prk_FC,IncreaseInGrossAssimilationRate_AP,"filled",'MarkerFaceAlpha','0.2','MarkerEdgeAlpha','0.2');
% Plot optimal point on top for comparison
hold on 
scatter3(change_aldolase_FC(2), change_prk_FC(2), IncreaseInGrossAssimilationRate_AP(2), 150,"pentagram",'filled','MarkerFaceColor',"#EDB120",'MarkerEdgeColor', 'black','MarkerFaceAlpha','1','MarkerEdgeAlpha','1');% plot optimal Aldolase on top in yellow
%scatter(change_prk_FC(2), IncreaseInGrossAssimilationRate_AP(2), 150,"pentagram",'filled','MarkerFaceColor',"#7E2F8E");% plot optimal PRK on top in purple
xticks([0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4]);
yticks([0 0.4 0.8 1.2 1.6 2.0 2.4 2.8]);
zticks([-0.4 -0.2 0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0]);
xlabel('Additional Aldolase','Rotation',10);
ylabel('Additional PRK','Rotation',-20);
zlabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
axis padded;

%% Scatter graph with a full set of 2000 nodes for SBPase-PRK pair
nexttile
% X is SBPase, Y is PRK, Z is Assimilation
scatter3(change_sbpase_FC,change_prk_FC,IncreaseInGrossAssimilationRate_SP,"filled",'MarkerFaceAlpha','0.2','MarkerEdgeAlpha','0.2');
% Plot optimal point on top for comparison
hold on 
scatter3(change_sbpase_FC(2), change_prk_FC(2), IncreaseInGrossAssimilationRate_SP(2), 150,"pentagram",'filled','MarkerFaceColor',"#EDB120",'MarkerEdgeColor', 'black','MarkerFaceAlpha','1','MarkerEdgeAlpha','1');% plot optimal SBPase top in green
%scatter(change_prk_FC(2), IncreaseInGrossAssimilationRate_SP(2), 150,"pentagram",'filled','MarkerFaceColor',"#7E2F8E");% plot optimal PRK on top in purple
xticks([0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4 4.8 5.2 5.6 6.0]);
yticks([0 0.4 0.8 1.2 1.6 2.0 2.4 2.8]);
zticks([-0.14 -0.04 0.06 0.16 0.26 0.36 0.46 0.56 0.66 0.76 0.86 0.96 1.06 1.16 1.26 1.36 1.46 1.56 1.66 1.76 1.86 1.96 2.06]);
xlabel('Additional SBPase','Rotation',10);
ylabel('Additional PRK','Rotation',-20);
zlabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
axis padded;
set(gcf, 'PaperOrientation', 'landscape');
print(gcf,fullfile('Outputs/rice_params/graphs',"sensitivity_two_enzymes_new.pdf"),'-dpdf','-fillpage');
%% Scatter graph with a full set of 2000 nodes for Rubisco-Aldolase-SBPase-PRK quad
%%First plotting only Rubisco against Assimilation
figure
q = tiledlayout(2,2,'Padding','compact');
title(q,'Four enzyme increases vs. changes in assimilation rate')
nexttile
scatter(change_rubisco_FC,IncreaseInGrossAssimilationRate_RSAP,"filled");
% Plot optimal point on top for comparison
hold on 
scatter(change_rubisco_FC(2), IncreaseInGrossAssimilationRate_RSAP(2), 150,"pentagram",'filled','MarkerFaceColor',"#D95319");% plot optimal Rubisco on top in red
xticks([0.0 0.2 0.4 0.6 0.8 1.0]);
yticks([-4.0 -2.0 0.0 2.0 4.0 6.0 8.0 10.0]);
xlabel('Additional Rubisco');
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
axis padded;

%%Then plotting only aldolase against assimilation
nexttile
scatter(change_aldolase_FC,IncreaseInGrossAssimilationRate_RSAP,"filled");
% Plot optimal point on top for comparison
hold on 
scatter(change_aldolase_FC(2), IncreaseInGrossAssimilationRate_RSAP(2), 150,"pentagram",'filled','MarkerFaceColor',"#EDB120");% plot optimal aldolase on top in yellow
xticks([0.0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4]);
yticks([-4.0 -2.0 0.0 2.0 4.0 6.0 8.0 10.0]);
xlabel('Additional Aldolase');
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
axis padded;

%%Then plotting only SBPase against assimilation
nexttile
scatter(change_sbpase_FC,IncreaseInGrossAssimilationRate_RSAP,"filled");
% Plot optimal point on top for comparison
hold on 
scatter(change_sbpase_FC(2), IncreaseInGrossAssimilationRate_RSAP(2), 150,"pentagram",'filled','MarkerFaceColor',"#77AC30");% plot optimal SBPase on top in green
xticks([0.0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4 4.8 5.2 5.6 6.0]);
yticks([-4.0 -2.0 0.0 2.0 4.0 6.0 8.0 10.0]);
xlabel('Additional SBPase');
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
axis padded;

%%Then plotting only PRK against assimilation
nexttile
scatter(change_prk_FC,IncreaseInGrossAssimilationRate_RSAP,"filled");
% Plot optimal point on top for comparison
hold on 
scatter(change_prk_FC(2), IncreaseInGrossAssimilationRate_RSAP(2), 150,"pentagram",'filled','MarkerFaceColor',"#7E2F8E");% plot optimal PRK on top in purple
xticks([0.0 0.4 0.8 1.2 1.6 2.0 2.4 2.8]);
yticks([-4.0 -2.0 0.0 2.0 4.0 6.0 8.0 10.0]);
xlabel('Additional PRK');
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
axis padded;
set(gcf, 'PaperOrientation', 'landscape');
print(gcf,fullfile('Outputs/rice_params/graphs',"sensitivity_four_enzymes_new.pdf"),'-dpdf','-fillpage');