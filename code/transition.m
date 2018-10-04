function updated_state_matrix = transition(def_path, prob_path, state_matrix, StateMatCols, state_mat_demog_group_idx, demog_group_def, DemogTblCols, param_names)
    
    dirlist = dir(strcat(def_path, '*.csv'));
    
    for file = dirlist'
            
        % Read in transition definition matrix
        [state_def_mat, StateDefCols] = read_table(strcat(def_path, file.name));

        % Find the people who are in the starting state within the
        % demographic group
        eligible_rows = state_mat_demog_group_idx;

        for field = fieldnames(StateDefCols)'
            field_string = string(field);
            start_state_value = state_def_mat(1, StateDefCols.(field_string));
            eligible_rows = find_indices(state_matrix, eligible_rows, StateMatCols.(field_string), '=', start_state_value); 
        end
                
        % Read in transition probabilities matrix
        [state_prob_mat, StateProbCols] = read_table(strcat(prob_path, file.name));
        
        % Get all param1 for the distributions
        param1_idx = find(contains(fieldnames(StateProbCols), param_names(1)));

        % Get all param2 for the distributions
        param2_idx = find(contains(fieldnames(StateProbCols), param_names(2)));

                    
        % Look up transition probability distribution parameters from table
        prob_row_idx = find_demog_rows(state_prob_mat, StateProbCols, demog_group_def, DemogTblCols, DemogTblCols);
        prob_params = state_prob_mat(prob_row_idx, :);
        
        % If a distribution is not provided, do not transition 
        if isempty(prob_params)
            continue
        end

        % Remove the starting state row (1st row) to make indexing
        % consistent with parameters
        state_def_mat = state_def_mat(2:end, :);
        
        % For each possible transition
        % this is in reference to the parameter index
        for i = 1:length(param1_idx)
            
            % Define the probability distribution for transitioning
            param1 = prob_params(param1_idx(i));
            param2 = prob_params(param2_idx(i));
            
            % Pick random probabilities from the demographic-specific
            % distribution
            transition_prob = betarnd(param1, param2, [1,length(eligible_rows)]);
            % Compare to random numbers picked from a uniform
            % distribution and determine who transitions
            unif_prob = rand([1, length(eligible_rows)]);
            to_transition = transition_prob > unif_prob;
                        
            state_mat_transition_rows = eligible_rows(to_transition);
            

            % Transition the people to the appropriate state
            for field = fieldnames(StateDefCols)'
                field_string = string(field);
                state_matrix(state_mat_transition_rows, StateMatCols.(field_string)) = state_def_mat(i, StateDefCols.(field_string));
            end
        end
    end
    updated_state_matrix = state_matrix;
end
