% Read in all 129 optimization results, average the results and export PR
% enzyme Vmax values to PR_constraints to serve as constraints for full optimizations

% Define the folder containing the Excel files
folder_129 = 'Outputs\Enzymes\129\';

% Get the list of text files containing enzyme Vmax of optimizations at 129 ppm
optimized_129_list = dir(fullfile(folder_129, '*.txt')); % get all text files in that folder

% Initialize an empty cell array to store the matrices
optimized_129_matrices = cell(numel(optimized_129_list), 1);

% Loop over the files
for i = 1:numel(optimized_129_list)
    % Read the current file using readmatrix
    optimized_129_matrix = readmatrix(fullfile(folder_129, optimized_129_list(i).name));
    
    % Store the matrix in the cell array
    optimized_129_matrices{i} = optimized_129_matrix;
end

% Combine all cells in the array 
combined_129_matrix = horzcat(optimized_129_matrices{:,1});
blank_rows = zeros(1,10);
combined_129_matrix_new = vertcat(combined_129_matrix(2:8,:), ...
    blank_rows,combined_129_matrix(9,:),blank_rows,combined_129_matrix(10:11,:),blank_rows,combined_129_matrix(12:25,:));

% Average across the columns for each enzyme Vmax
avg_optimized_129 = mean(combined_129_matrix,2);

% Extract rows for photorespiratory constraints and save in text file
PR_constraints = avg_optimized_129(12:18);
writematrix(PR_constraints,'PR_constraints.txt');

% Update results in file
% writematrix(combined_129_matrix_new,'Outputs/rice_params/Results_optimization_rice_new_2.xlsx','Sheet','129','Range','K4')
