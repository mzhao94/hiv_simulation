% What are the dimensions I need in order to measure prevalence? 
% 1) I need the variable I want the prevalence of
% 2) I need the dimension I want to take the prevalence over
    % for now we're going to only do this over time
% 3) I need the population characteristics we're trying to calculate the
% metric over
    % The ability to graph multiple groups over time should be built in

% Need to build in functionality for population_def indicating that there
% should be no restriction on population
function results = prevalence_by_time(state_matrices_path, StateMatCols, main_var, main_var_values, population_def, PopMatCols)
    
    % Empty list to store results in
    results = [];
    
    % Get all file names in the state_matrices_path
    state_matrix_files = dir(strcat(state_matrices_path, '*.mat'));
    state_matrix_names = natsortfiles({state_matrix_files.name});
    
    for file = state_matrix_names
        state_matrix_struct = load(strcat(state_matrices_path, string(file)));
        state_matrix = state_matrix_struct.state_matrix;
        
        if size(population_def, 1) > 0
            disp(population_def)
            disp(size(population_def, 1))
            % Find the people we are calculating prevalence over
            demog_mat_indices = find_demog_rows(state_matrix, StateMatCols, population_def, PopMatCols, PopMatCols);
            state_matrix = state_matrix(demog_mat_indices, main_var);
        end
        % Calculate prevalence of the main_var in the subsetted
        % state_matrix
        eligible_indices = find_indices(state_matrix, 1:size(state_matrix, 1), StateMatCols.(main_var), '=', main_var_values);
        count = size(eligible_indices, 2);
        disp(count)
        prevalence = count/size(state_matrix,1);
        results = [results, prevalence];
    end  
 