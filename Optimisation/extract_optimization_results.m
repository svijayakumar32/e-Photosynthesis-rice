% Read in optimization results for each Cc (130-380) and average the results one Cc at a time
% Define the folder containing the results
folder_130 = 'Outputs\Enzymes\130\';

% Get the list of text files containing enzyme Vmax of optimizations at 130 ppm
optimized_130_list = dir(fullfile(folder_130, '*.txt')); % get all text files in that folder

% Initialize an empty cell array to store the matrices
optimized_130_matrices = cell(numel(optimized_130_list), 1);

% Loop over the files
for i = 1:numel(optimized_130_list)
    % Read the current file using readmatrix
    optimized_130_matrix = readmatrix(fullfile(folder_130, optimized_130_list(i).name));
    
    % Store the matrix in the cell array
    optimized_130_matrices{i} = optimized_130_matrix;
end

% Combine all cells in the array 
combined_130_matrix = horzcat(optimized_130_matrices{:,1});
blank_rows = zeros(1,10);
combined_130_matrix_new = vertcat(combined_130_matrix(2:8,:), ...
    blank_rows,combined_130_matrix(9,:),blank_rows,combined_130_matrix(10:11,:),blank_rows,combined_130_matrix(12:25,:));
% Average across the columns for each enzyme Vmax
avg_optimized_130 = mean(combined_130_matrix_new,2);

% for j = 2:numel(combined_130_matrix_new)
%     % Calculate the protein amount
%     combined_130_protein{j} = (combined_130_matrix_new(j)/MWKcat(j,2))*(MWKcat(j,3)/1000);
% end
