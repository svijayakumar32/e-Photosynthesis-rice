%% Plot A/Cc fits using limitations file
%% Load data
% Load limitations
lim = readtable("completeFit_rice_limitations.csv"); % Limitation data
ccVars = lim.Properties.VariableNames(startsWith(lim.Properties.VariableNames, 'Cc')); % Cc columns
aVars  = lim.Properties.VariableNames(startsWith(lim.Properties.VariableNames, 'A')); % A columns

% Load Gstar file
gstar = readtable("Gstars_Pa_rice.csv");

% Extract Gstar values
Gstars = zeros(8,1);
for i=1:numel(Gstars)
    Gstars(i)= gstar{i,2};
end

% Extract data points by reading data from tables
data = cell(1, 8);
for j = 1:8
    filename = sprintf('completeFits_rice_dataTable_%d.csv', j);
    data{j} = readtable(filename);
end

%% Graph plotting - group of columns in lim
nGroups = min(numel(ccVars), numel(aVars)) / 3;  % Columns per group
nPlots = 3;     % Number of sets to plot (Cc/Ac, Cc/Aj, Cc/Ap)

% Specify 4x2 layout for panel figure 
fig = figure;
t = tiledlayout(fig,4, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

for k = 0:(nGroups-1)
    ax = nexttile;
    
    % Column names for each group k
    ccCols = ccVars(k*3 + (1:3));
    aCols  = aVars(k*3 + (1:3));

    hold(ax, 'on');

    % Plot limitation sets using column indices
    p = scatter(ax,lim.(ccCols{1}), lim.(aCols{1}), 'filled', 'MarkerFaceColor', "#C00000");
    q = scatter(ax,lim.(ccCols{2}), lim.(aCols{2}), 'filled', 'MarkerFaceColor', "#0000FF");
    r = scatter(ax,lim.(ccCols{3}), lim.(aCols{3}), 'filled', 'MarkerFaceColor', "#FFC000");

    % Plot data using column indices
    d = scatter(ax, data{1, k+1}.Cc, data{1, k+1}.A, 'filled', 'MarkerFaceColor', 'black');

    % Plot Gstar line 
    g = plot(ax, [Gstars(k+1) Gstars(k+1)], [-20 60], '--', 'Color', '#7030A0');

    % Set axis limits, labels and legend
    xlim(ax, [0 250]);
    ylim(ax, [-20 60]);
    ax.FontSize = 10;
    xlabel(ax, 'C_{c} (Pa)', 'FontSize', 12);
    ylabel(ax, 'A (μmol m^{−2} s^{−1})', 'FontSize', 12);

    legend(ax, [p, q, r, g, d], {'\it{A_{c}}', '\it{A_{j}}', '\it{A_{p}}', '\Gamma^{*}', 'Data'}, ...
           'Location', 'southeast', 'FontSize', 10);
    
    hold(ax, 'off');
end

% Save figure in pdf file
print(fig, 'Curve_fitting_A_Cc_new.pdf', '-dpdf', '-fillpage');
