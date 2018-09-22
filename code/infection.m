% This function will take in a particular demographic group, a state
% matrix, and a mixing matrix and will transition people in this
% demographic group to infection based on a uniform random number draw

function updated_state_matrix = infection(reference_demog_group_def, demog_table, DemogTblCols, mixing_matrix, state_matrix, StateMatCols)
    
    % Read in the betas for this particular demographic group (i.e. row in
    % the matrix)
    betas = get_mixing_matrix_row(mixing_matrix, reference_demog_group_def, DemogTblCols);
    
    % Placeholder for adding up all the infection probabilities over all
    % demographic group our reference group is mixing with
    total_infec_prob = 0;
    
    % ASSUMPTION: the demog table and the mixing matrix have
    % the same order of demographic groups
    
    % For each of all possible demographic groups
    for demog_group_idx = 1:size(betas)
        
        % Get the demographic group definition we are mixing with
        other_demog_group_def = demog_table(demog_group_idx,:);
        
        % Calculate infection probability
        beta = betas(demog_group_idx);
        infec_prob = beta * calc_infec_prob(state_matrix, StateMatCols, other_demog_group_def, DemogTblCols);
        
        % Add this to the total infection probability over all demographic
        % groups we are mixing with
        total_infec_prob = total_infec_prob + infec_prob;
        
    
    % Find the people in the demographic group
    reference_demog_mat_indices = find_demog_rows(state_matrix, StateMatCols, reference_demog_group_def, demog_col_struct);
    
    % Find eligible people to infect in our reference demographic group (i.e. no hiv, no prep, no art, alive)
    eligible_rows = find_indices(state_matrix, reference_demog_mat_indices, 'alive', '=', 1);
    eligible_rows = find_indices(state_matrix, eligible_rows, 'hiv', '=', 0);
    eligible_rows = find_indices(state_matrix, eligible_rows, 'prep', '=', 0);
    
    % Do a random number draw for all of the eligible rows
    unif_prob = rand([1, length(eligible_rows)]);

    % Compare to random numbers picked from a uniform
    % distribution and determine who is infected
    to_infect = total_infec_prob > unif_prob;
                        
    % Get indices to infect based on people who are eligible
    state_mat_infec_rows = eligible_rows(to_infect);
    
    % Infect chosen people
    state_matrix(state_mat_infec_rows, StateMatCols.hiv) = 1;
    
    % Return the updated state matrix
    updated_state_matrix = state_matrix;
    
    end
end