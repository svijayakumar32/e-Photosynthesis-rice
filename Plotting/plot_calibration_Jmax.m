%% Plot and examine results of the parameter calibration to identify the best scaling factors αEnzymes and αRubisco

%% Load Jmax data
load('Jmax_simple_3_result.mat');

% Create figure
figure;
hold on;

% Specify Cc and Gross Assimilation data for FvCB model
FvCB_x = Farq_Matrix(:, 5);
FvCB_y = Farq_Matrix(:, 4);
% Plot FvCB points with markers and dashed lines
plot(FvCB_x, FvCB_y, 'o', ...   % '--o' specifies dashed line with circle markers
     'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8);  % black markers

% For e_Photosynthesis data, we are plotting multiple lines, each associated with a scaling factor between 0.52-1.5 (n=50)
% Each scaling factor has 35 Cc values (data points) associated with Jmax 
chunkSize = size(numbers,1);

% Generate different colors for each new scaling factor
numChunks = ceil(size(ePhoto_Matrix, 1) / chunkSize);
colors = lines(numChunks);

% Convert scaling factors into string array for a data table/legend, including FvCB as first series
labels = string(SSR_Matrix(:, 1));
new_labels = vertcat("FvCB", labels);
Cc_labels = string(Farq_Matrix(1:chunkSize,5));
ePhoto_table = array2table(numbers, 'VariableNames', labels,RowNames = Cc_labels);

% If choosing to only plot a subset of scaling factors 
% scaling_idx = [4:5:49]; % 0.58-1.48

% Plot Cc by Gross Assimilation for ePhotosynthesis model
for chunkIdx = 1:numChunks
    % Calculate row indices for the current chunk
    startRow = (chunkIdx - 1) * chunkSize + 1;
    endRow = min(chunkIdx * chunkSize, size(ePhoto_Matrix, 1));  

    % Extract x and y data for the current chunk
    chunk_x = ePhoto_Matrix(startRow:endRow, 5);  % Cc (x)
    chunk_y = ePhoto_Matrix(startRow:endRow, 4);  % Gross Assimilation (y)
    
    % Plot Cc (x) against Gross Assimilation (y) with specified colours for each scaling factor as dashed lines
    plot(chunk_x, chunk_y, '--o', 'Color', colors(chunkIdx, :));

    % Plot each row in the current chunk
    % for rowIdx = startRow:endRow
    %     % Plot Cc (x) against Gross Assimilation (y) with specified colours for each scaling factor
    %     plot(ePhoto_Matrix(rowIdx, 5), ePhoto_Matrix(rowIdx, 4),'--o','Color', colors(chunkIdx, :));
    % end
end

hold off;

% Set the x-axis limits from 180 to 640
xlim([180 640]);

% Add labels, title, and legend
xlabel('C_c');
ylabel('Gross Assimilation');
title('Scaling Factors for RuBP Regeneration');

% Append labels as the legend
legend(new_labels, 'Location', 'eastoutside');

% Print in pdf
set(gcf, 'PaperOrientation', 'landscape');
print(gcf,fullfile('Outputs/rice_params/graphs',"Scaling_Factors_for_RuBP_Regeneration"),'-dpdf','-fillpage');
