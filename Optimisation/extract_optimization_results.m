% Read in optimization results for each Cc (130-380) and average the results - can I
% automate this for all Cc levels at one?

% One Cc at a time...
% Define the folder containing the Excel files
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

% Update results in file
writematrix(combined_130_matrix_new,'Results_optimization_rice_new_2.xlsx','Sheet','130','Range','I4')

% %% Convert Vmax to protein (WORK IN PROGRESS)
% % Load kcat and MW data
% Enzyme=importdata('MW&Kcat.txt');
% MWKcat=Enzyme.data;
% %MWKcat([7,9,12],:) = []; %remove rows corresponding to V8, V10 and V16
% % Create empty cell to store protein concentrations
% combined_130_matrix_new(1,:) = [];
% combined_130_protein = cell(numel(optimized_130_list), 1);
% % Loop over the rows to obtain
%     combined_130_protein{1} = (combined_130_matrix_new(1)/MWKcat(1,2))*(MWKcat(1,3)/1000)*1.1;
% for j = 2:numel(combined_130_matrix_new)
%     % Calculate the protein amount
%     combined_130_protein{j} = (combined_130_matrix_new(j)/MWKcat(j,2))*(MWKcat(j,3)/1000);
% end
