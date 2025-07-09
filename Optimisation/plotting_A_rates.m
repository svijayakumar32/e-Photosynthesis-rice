% Plot assimilation rates for rice data

% Import assimilation rates from results file
file = 'Outputs/rice_params/Results_optimization_rice_new_3.xlsx';
optimized_average_As = 'G4'; % average assimilation rate from each sheet
sheet_names = (string(num2cell([130:10:380])));

% Import non-optimized and scenario-based assimilation rates from text file (generated after running CalculateGrossAssimilation)
non_optimized_A = num2cell(readmatrix('non_optimized_A.txt'));
reoptimized_low = num2cell(readmatrix('low_A.txt'));
reoptimized_ambient = num2cell(readmatrix('ambient_A.txt'));
reoptimized_elevated = num2cell(readmatrix('elevated_A.txt'));
twofold_low = num2cell(readmatrix('low2_A.txt'));
twofold_ambient = num2cell(readmatrix('ambient2_A.txt'));
twofold_elevated = num2cell(readmatrix('elevated2_A.txt'));

% Create cell array to import average assimilation rates from each Ci/sheet
reoptimized_assimilation_data = cell(8, numel(sheet_names));%non opt, opt and three strategies (2/4/6 optimized enzymes) and three strategies with twofold increases in Vmax
%reoptimized_assimilation_data = cell(2, numel(sheet_names)); % only non opt and opt A

% Convert sheet_names into x values 
Cc = str2double(sheet_names);

% Loop through all sheets and import values into assimilation_data
for i = 1:numel(Cc)
    reoptimized_assimilation_data{1,i} = non_optimized_A{i,1};
    reoptimized_assimilation_data{2,i} = xlsread(file, char(sheet_names(i)), optimized_average_As);
    reoptimized_assimilation_data{3,i} = reoptimized_low{i,1};
    reoptimized_assimilation_data{4,i} = reoptimized_ambient{i,1};
    reoptimized_assimilation_data{5,i} = reoptimized_elevated{i,1};
	reoptimized_assimilation_data{6,i} = twofold_low{i,1};
    reoptimized_assimilation_data{7,i} = twofold_ambient{i,1};
    reoptimized_assimilation_data{8,i} = twofold_elevated{i,1};
end
    
% Percentage change in A
percent_change_optimized=vertcat(sheet_names,(cell2mat(reoptimized_assimilation_data(2,:))-cell2mat(reoptimized_assimilation_data(1,:)))./cell2mat(reoptimized_assimilation_data(1,:))*100);
percent_change_low=vertcat(sheet_names,(cell2mat(reoptimized_assimilation_data(3,:))-cell2mat(reoptimized_assimilation_data(1,:)))./cell2mat(reoptimized_assimilation_data(1,:))*100);
percent_change_ambient=vertcat(sheet_names,(cell2mat(reoptimized_assimilation_data(4,:))-cell2mat(reoptimized_assimilation_data(1,:)))./cell2mat(reoptimized_assimilation_data(1,:))*100);
percent_change_elevated=vertcat(sheet_names,(cell2mat(reoptimized_assimilation_data(5,:))-cell2mat(reoptimized_assimilation_data(1,:)))./cell2mat(reoptimized_assimilation_data(1,:))*100);
percent_change_low2=vertcat(sheet_names,(cell2mat(reoptimized_assimilation_data(6,:))-cell2mat(reoptimized_assimilation_data(1,:)))./cell2mat(reoptimized_assimilation_data(1,:))*100);
percent_change_ambient2=vertcat(sheet_names,(cell2mat(reoptimized_assimilation_data(7,:))-cell2mat(reoptimized_assimilation_data(1,:)))./cell2mat(reoptimized_assimilation_data(1,:))*100);
percent_change_elevated2=vertcat(sheet_names,(cell2mat(reoptimized_assimilation_data(8,:))-cell2mat(reoptimized_assimilation_data(1,:)))./cell2mat(reoptimized_assimilation_data(1,:))*100);
%% Plot assimilation rates' averages against Cc for strategies using optimized Vmax
figure;
% Plot scatter graph of Cc (x) against average assimilation rates
% Plot non-optimized A
p1 = scatter(Cc, cell2mat(reoptimized_assimilation_data(1,:)),100,'black','filled','o','MarkerEdgeColor','black','LineWidth',0.5); 
p1.DisplayName = 'Non-Optimised';
hold on
% Plot optimized A averages
p2 = scatter(Cc, cell2mat(reoptimized_assimilation_data(2,:)),100,'white','filled','^','MarkerEdgeColor','black','LineWidth',0.5);
p2.DisplayName = 'Optimised';
hold on
% Plot non-optimized A with optimized A averages for Aldolase, SBPase, PRK at 250
p4 = scatter(Cc, cell2mat(reoptimized_assimilation_data(4,:)),100,'magenta','filled','square');
p4.DisplayName = 'Current';
hold on
% Plot non-optimized A with optimized A averages for Aldolase, SBPase, PRK, FBPase, PGK at 360
p5 = scatter(Cc, cell2mat(reoptimized_assimilation_data(5,:)),100,'blue','filled','p');
p5.DisplayName = 'Future';
hold on
% Plot non-optimized A with optimized A averages for SBPase at 130
p3 = scatter(Cc, cell2mat(reoptimized_assimilation_data(3,:)),100,'red','filled','h');
p3.DisplayName = 'Stress';
hold on
xlim([130,380]);
xticks([130:50:380]);
xticklabels([130:50:380]);
ylim([0 50]);
yticks([0 10 20 30 40 50]);
xlabel('C_{c} (μmol mol^{-1})');
ylabel({'Gross CO_2 Assimilation', '(μmol m^{-2} s^{-1})'});
ax = gca;
ax.FontSize = 20;
legend([p1, p2, p3, p4, p5],'Location', 'southeast','FontSize', 12)
hold off
set(gcf, 'PaperOrientation', 'landscape');
%print(gcf,fullfile('Outputs/rice_params/graphs',"Cc_vs_Assimilation_new_strategies"),'-dpdf','-bestfit');
%exportgraphics(gcf, fullfile('Outputs/rice_params/graphs', "Cc_vs_Assimilation_new_strategies3.pdf"), 'ContentType', 'vector', 'Resolution', 300, 'ColorSpace', 'rgb');

%% Plot assimilation rates' averages against Cc for strategies using twofold increases in Vmax
figure;
% Plot scatter graph of Cc (x) against average assimilation rates
% Plot non-optimized A
scatter(Cc, cell2mat(reoptimized_assimilation_data(1,:)),100,'black','filled','o','MarkerEdgeColor','black','LineWidth',0.5); 
hold on
% Plot optimized A averages
scatter(Cc, cell2mat(reoptimized_assimilation_data(2,:)),100,'white','filled','^','MarkerEdgeColor','black','LineWidth',0.5);
hold on
% Plot non-optimized A with twofold increases in SBPase, Aldolase, PRK
scatter(Cc, cell2mat(reoptimized_assimilation_data(7,:)),100,'magenta','filled','square');
hold on
% Plot non-optimized A with twofold increases in SBPase, Aldolase, PRK, FBPase, TK
scatter(Cc, cell2mat(reoptimized_assimilation_data(8,:)),100,'blue','filled','p');
hold on
% Plot non-optimized A with twofold increases in SBPase
scatter(Cc, cell2mat(reoptimized_assimilation_data(6,:)),100,'red','filled','h');
hold on
xlim([130,380]);
xticks([130:50:380]);
xticklabels([130:50:380]);
ylim([0 50]);
yticks([0 10 20 30 40 50]);
xlabel('C_{c} (μmol mol^{-1})');
ylabel({'Gross CO_{2} Assimilation', '(μmol m^{-2} s^{-1})'});
ax = gca;
ax.FontSize = 20;
hold off
set(gcf, 'PaperOrientation', 'landscape');
%print(gcf,fullfile('Outputs/rice_params/graphs',"Cc_vs_Assimilation_new_strategies_twofold"),'-dpdf','-bestfit');
%exportgraphics(gcf, fullfile('Outputs/rice_params/graphs', "Cc_vs_Assimilation_new_strategies_twofold3.pdf"), 'ContentType', 'vector', 'Resolution', 300, 'ColorSpace', 'rgb');

