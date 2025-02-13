%% Load data
load output_enzyme_adjustment_test_250.mat; 
%% Create outputs for increases in gross CO2 assimilation rate
[IncreaseInGrossAssimilationRate_RSAP,... %four
    IncreaseInGrossAssimilationRate_RA,IncreaseInGrossAssimilationRate_RS,IncreaseInGrossAssimilationRate_RP,...%two
        IncreaseInGrossAssimilationRate_AS,IncreaseInGrossAssimilationRate_AP,IncreaseInGrossAssimilationRate_SP] = deal(zeros(2000,1));
%% Convert values for plotting
%Expressing the y-axis as the increases/decreases in net CO2 assimilation rate
for n=1:length(IncreaseInGrossAssimilationRate_RSAP)
IncreaseInGrossAssimilationRate_RSAP(n) = GrossAssimilationRate_RSAP(n)-GrossAssimilationRate_RSAP(1);
IncreaseInGrossAssimilationRate_RA(n) = GrossAssimilationRate_RA(n)-GrossAssimilationRate_RA(1);
IncreaseInGrossAssimilationRate_RS(n) = GrossAssimilationRate_RS(n)-GrossAssimilationRate_RS(1);
IncreaseInGrossAssimilationRate_RP(n) = GrossAssimilationRate_RP(n)-GrossAssimilationRate_RP(1);
IncreaseInGrossAssimilationRate_AS(n) = GrossAssimilationRate_AS(n)-GrossAssimilationRate_AS(1);
IncreaseInGrossAssimilationRate_AP(n) = GrossAssimilationRate_AP(n)-GrossAssimilationRate_AP(1);
IncreaseInGrossAssimilationRate_SP(n) = GrossAssimilationRate_SP(n)-GrossAssimilationRate_SP(1);
end

%Expressing the x-axis as the multiplicative change in protein content
change_aldolase_FC = aldolase_FC-1;
change_rubisco_FC = rubisco_FC-1;
change_sbpase_FC = sbpase_FC-1;
change_prk_FC = prk_FC-1;
%% Create scatter plots for each enzyme pair
%% Rubisco
figure;
%scatter(rubisco_FC, GrossAssimilationRate_R, 'filled','MarkerFaceColor',"#D95319");
scatter(change_rubisco_FC, IncreaseInGrossAssimilationRate_RA, 'filled','MarkerFaceColor',"#D95319");
scatter(change_rubisco_FC, IncreaseInGrossAssimilationRate_RS, 'filled','MarkerFaceColor',"#D95319");
scatter(change_rubisco_FC, IncreaseInGrossAssimilationRate_RP, 'filled','MarkerFaceColor',"#D95319");

% Plot optimal point for comparison
hold on 
%scatter(rubisco_FC(1), GrossAssimilationRate_R(1), 75,"v", 'filled','MarkerFaceColor',"#0072BD");% plot min FC (=1) on top
scatter(change_rubisco_FC(3), IncreaseInGrossAssimilationRate_R(3), 150,"pentagram",'filled','MarkerFaceColor',"#0072BD");% plot optimal point on top
%scatter(rubisco_FC(2), GrossAssimilationRate_R(2), 75,"^", 'filled','MarkerFaceColor',"#0072BD");% plot max FC on top
% Add axes labels
%xlabel('Rubisco Fold Change');
%ylabel('Gross Assimilation Rate (μmol m^{−2} s^{−1})');
xlabel('Additional Rubisco');
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
xlim([0 1.0]);
ylim([-4.5 0]);
xticks([0 0.2 0.4 0.6 0.8 1.0]);
yticks([-4.5 -4.0 -3.5 -3.0 -2.5 -2.0 -1.5 -1.0 -0.5 0]);
legend('','optimalFC','Location', 'northeastoutside');

grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
% Save the plot as an image file
saveas(gcf, 'scatter_R_new.png');
saveas(gcf, 'scatter_R_new.pdf');
%% Aldolase
figure;
scatter(change_aldolase_FC, IncreaseInGrossAssimilationRate_A, 'filled','MarkerFaceColor',"#0072BD");
% Plot optimal point for comparison
hold on 
%scatter(aldolase_FC(1), GrossAssimilationRate_A(1), 75,"v", 'filled','MarkerFaceColor',"#EDB120"); %min
scatter(change_aldolase_FC(3), IncreaseInGrossAssimilationRate_A(3), 150,"pentagram", 'filled','MarkerFaceColor',"#EDB120"); %opt
%scatter(aldolase_FC(2), GrossAssimilationRate_A(2), 75,"^", 'filled','MarkerFaceColor',"#EDB120"); %max

% Add axes labels
%xlabel('Aldolase Fold Change');
%ylabel('Gross Assimilation Rate (μmol m^{−2} s^{−1})');
xlabel('Additional Aldolase');
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
xlim([0 2.0]);
ylim([0 0.25]);
xticks([0 0.4 0.8 1.2 1.6 2.0]);
yticks([0 0.05 0.10 0.15 0.20 0.25]);
legend('','optimalFC','Location', 'northeastoutside');

grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
% Save the plot as an image file
saveas(gcf, 'scatter_A_new.png');
saveas(gcf, 'scatter_A_new.pdf');

%% SBPase
figure;
scatter(change_sbpase_FC, IncreaseInGrossAssimilationRate_S, 'filled','MarkerFaceColor',"#EDB120");
% Plot optimal point for comparison
hold on 
%scatter(sbpase_FC(1), GrossAssimilationRate_S(1),75,"v",'filled','MarkerFaceColor',"#0072BD");%min
scatter(change_sbpase_FC(3), IncreaseInGrossAssimilationRate_S(3),150,"pentagram",'filled','MarkerFaceColor',"#0072BD");%opt
%scatter(sbpase_FC(2), GrossAssimilationRate_S(2),75,"^",'filled','MarkerFaceColor',"#0072BD");%max

% Add axes labels
%xlabel('SBPase Fold Change');
%ylabel('Gross Assimilation Rate (μmol m^{−2} s^{−1})');
xlabel('Additional SBPase');
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
xlim([0 4.8]);
ylim([0 0.8]);
xticks([0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4 4.8]);
yticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8]);
legend('','optimalFC','Location', 'northeastoutside');

grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
% Save the plot as an image file
saveas(gcf, 'scatter_S_new.png');
saveas(gcf, 'scatter_S_new.pdf');

%% Rubisco vs Gross Assimilation RSAP
figure;
scatter(change_rubisco_FC, IncreaseInGrossAssimilationRate_RSA, 'filled','MarkerFaceColor',"#D95319");
% Plot optimal point for comparison
hold on
%scatter(rubisco_FC(1), GrossAssimilationRate_RSA(1), 75,"v", 'filled','MarkerFaceColor',"#0072BD");% plot min FC (=1) on top
scatter(change_rubisco_FC(3), IncreaseInGrossAssimilationRate_RSA(3), 150,"pentagram",'filled','MarkerFaceColor',"#0072BD");% plot optimal point on top
%scatter(rubisco_FC(2), GrossAssimilationRate_RSA(2), 75,"^", 'filled','MarkerFaceColor',"#0072BD");% plot max FC on top
% Add axes labels
%xlabel('Rubisco Fold Change');
%ylabel('Gross Assimilation Rate (μmol m^{−2} s^{−1})');
xlabel('Additional Rubisco');
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
xlim([0 1.0]);
ylim([0 9.0]);
xticks([0 0.2 0.4 0.6 0.8 1.0]);
yticks([0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0]);
legend('','optimalFC','Location', 'northeastoutside');

grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
% Save the plot as an image file
saveas(gcf, 'scatter_R_vs_RSA_new.png');
saveas(gcf, 'scatter_R_vs_RSA_new.pdf');

%% Aldolase vs Gross Assimilation RSA
figure;
scatter(change_aldolase_FC, IncreaseInGrossAssimilationRate_RSA, 'filled','MarkerFaceColor',"#0072BD");
% Plot optimal point for comparison
hold on 
%scatter(aldolase_FC(1), GrossAssimilationRate_RSA(1), 75,"v", 'filled','MarkerFaceColor',"#EDB120"); %min
scatter(change_aldolase_FC(3), IncreaseInGrossAssimilationRate_RSA(3), 150,"pentagram", 'filled','MarkerFaceColor',"#EDB120"); %opt
%scatter(aldolase_FC(2), GrossAssimilationRate_RSA(2), 75,"^", 'filled','MarkerFaceColor',"#EDB120"); %max

% Add axes labels
%xlabel('Aldolase Fold Change');
%ylabel('Gross Assimilation Rate (μmol m^{−2} s^{−1})');
xlabel('Additional Aldolase');
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
xlim([0 2.0]);
ylim([0 9.0]);
xticks([0 0.4 0.8 1.2 1.6 2.0]);
yticks([0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0]);
legend('','optimalFC','Location', 'northeastoutside');

grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
% Save the plot as an image file
saveas(gcf, 'scatter_A_vs_RSA_new.png');
saveas(gcf, 'scatter_A_vs_RSA_new.pdf');

%% SBPase vs Gross Assimilation RSA
figure;
scatter(change_sbpase_FC, IncreaseInGrossAssimilationRate_RSA, 'filled','MarkerFaceColor',"#EDB120");
% Plot optimal point for comparison
hold on 
%scatter(sbpase_FC(1), GrossAssimilationRate_RSA(1),75,"v",'filled','MarkerFaceColor',"#0072BD");%min
scatter(change_sbpase_FC(3), IncreaseInGrossAssimilationRate_RSA(3),150,"pentagram",'filled','MarkerFaceColor',"#0072BD");%opt
%scatter(sbpase_FC(2), GrossAssimilationRate_RSA(2),75,"^",'filled','MarkerFaceColor',"#0072BD");%max

% Add axes labels
%xlabel('SBPase Fold Change');
%ylabel('Gross Assimilation Rate (μmol m^{−2} s^{−1})');
xlabel('Additional SBPase');
ylabel('ΔCO_{2} uptake (μmol m^{−2} s^{−1})');
xlim([0 4.8]);
ylim([0 9.0]);
xticks([0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4 4.8]);
yticks([0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0]);
legend('','optimalFC','Location', 'northeastoutside');

grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data
% Save the plot as an image file
saveas(gcf, 'scatter_S_vs_RSA_new.png');
saveas(gcf, 'scatter_S_vs_RSA_new.pdf');

%% Rubisco, Aldolase and SBPase
% Create a 3D scatter plot for RSA result
figure;
scatter3(change_rubisco_FC, change_aldolase_FC, change_sbpase_FC, [], IncreaseInGrossAssimilationRate_RSA, 'filled'); % 'filled' adds color to the markers based on 'y'
colorbar; % Display a colorbar to show the mapping of 'y' values to colors

% Add axes labels
xlabel('Additional Rubisco');
ylabel('Additional Aldolase');
zlabel('Additional SBPase');
%title('3D Scatter Plot');

% Customize a colormap if needed
colormap('winter'); % You can change the colormap (e.g., 'jet', 'parula', 'viridis', etc.)
legend({'ΔCO_{2} uptake (μmol m^{−2} s^{−1})'}, 'Location', 'northeastoutside'); % 'Location' can be adjusted as needed

% Adjust the view angle if needed
view(45, 30); % Example view angle (azimuth, elevation)

% Show the plot
grid on; % Add grid lines
axis tight; % Fit the axes tightly to the data

% Save the plot as an image file
saveas(gcf, 'scatter_RSA_new.png');
saveas(gcf, 'scatter_RSA.pdf');
