% This function will take in a particular demographic group, a state
% matrix, and a mixing matrix and will transition people in this
% demographic group to infection based on a uniform random number draw

function updated_state_matrix = infection(reference_demog_group_def, DemogTblCols, mixing_matrix, mixing_table, MixingTblCols, state_matrix, StateMatCols)
    
    % Returns the row index based on the mixing_mat_def
    % demographic group table
    demog_row = find_demog_rows(mixing_table, MixingTblCols, reference_demog_group_def, DemogTblCols, MixingTblCols);

    % Based on the demographic row, extract the betas from the
    % corresponding row number in the mixing matrix
    betas = mixing_matrix(demog_row,:);
        
    % Placeholder for adding up all the infection probabilities over all
    % demographic group our reference group is mixing with
    total_infec_prob = 0;
    
    % For each of all possible demographic groups
    for demog_group_idx = 1:size(betas, 2)

        % Get the demographic group definition we are mixing with
        other_demog_group_def = mixing_table(demog_group_idx,:);

        % Calculate infection probability
        beta = betas(demog_group_idx);
        infec_prob = beta * calc_infec_prob(state_matrix, StateMatCols, other_demog_group_def, DemogTblCols, MixingTblCols);
        
        % Add this to the total infection probability over all demographic
        % groups we are mixing with
        total_infec_prob = total_infec_prob + infec_prob;
    
    end  

    % Find the people in the demographic group
    reference_demog_mat_indices = find_demog_rows(state_matrix, StateMatCols, reference_demog_group_def, DemogTblCols, MixingTblCols);
    
    % Find eligible people to infect in our reference demographic group (i.e. no hiv, no prep, no art, alive)
    susceptible = find_indices(state_matrix, reference_demog_mat_indices, StateMatCols.alive, '=', 1);
    susceptible = find_indices(state_matrix, susceptible, StateMatCols.hiv, '=', 0);
    susceptible = find_indices(state_matrix, susceptible, StateMatCols.prep, '=', 0);
    susceptible = find_indices(state_matrix, susceptible, StateMatCols.art, '=', 0);

   
    % (I/N) * S
    total_infec_prob = total_infec_prob * length(susceptible);
    
    % Do a random number draw for all of the eligible rows
    unif_prob = rand([1, length(susceptible)]);

    % Compare to random numbers picked from a uniform
    % distribution and determine who is infected
    to_infect = total_infec_prob > unif_prob;
                        
    % Get indices to infect based on people who are eligible
    state_mat_infec_rows = susceptible(to_infect);
    
    % Infect chosen people
    state_matrix(state_mat_infec_rows, StateMatCols.hiv) = 1;
    
    % Return the updated state matrix
    updated_state_matrix = state_matrix;
    
end