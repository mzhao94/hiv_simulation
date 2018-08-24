% Run the simulation parameters
hiv_simulation_parameters;

% Read in the initial population from a csv
% Assumption: the csv file has all the people pre-populated with
% characteristics
[state_matrix, StateMatCols] = read_table(init_pop_file);

% Read in all possible demographic groups
[demog_table, DemogTblCols] = read_table(demog_groups_file);

% Create array for case when everyone is eligible for subsetting
all_eligible = 1:size(state_matrix, 1);


% For each time period
for t = 1:T
    
    % For each demographic group
    for demog_group_def = demog_table.'
        
        % Get row indices for people in the state matrix in this
        % demographic group
        state_mat_demog_group_idx = find_demog_rows(state_matrix, StateMatCols, demog_group_def, DemogTblCols);

        
        % If there's no one in the state matrix in this demographic group, continue to the next
        % one
        if isempty(state_mat_demog_group_idx)
            continue
        end 
                        
        %%%%%%%%%%%%%%%%%%%%%%%%%% BIRTHS %%%%%%%%%%%%%%%%%%%%%%%%%%
        state_matrix = transition(birth_transition_definition_path, birth_transition_probabilities_path, ...
                                    state_matrix, StateMatCols, state_mat_demog_group_idx, ...
                                    demog_group_def, DemogTblCols, param_names);
                       

        %%%%%%%%%%%%%%%%%%%%%%%%%% DEATHS %%%%%%%%%%%%%%%%%%%%%%%%%% 
        state_matrix = transition(death_transition_definition_path, death_transition_probabilities_path, ...
                                    state_matrix, StateMatCols, state_mat_demog_group_idx, ...
                                    demog_group_def, DemogTblCols, param_names);
        

                                
        %%%%%%%%%%%%%%%%%%%%%%%%%% HIV STATUS TRANSITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%
        state_matrix = transition(hiv_transition_definition_path, hiv_transition_probabilities_path, ...
                                    state_matrix, StateMatCols, state_mat_demog_group_idx, ...
                                    demog_group_def, DemogTblCols, param_names);
        
                                
        %%%%%%%%%%%%%%%%%%%%%%%%%% ART STATUS TRANSITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%
        state_matrix = transition(art_transition_definition_path, art_transition_probabilities_path, ...
                            state_matrix, StateMatCols, state_mat_demog_group_idx, ...
                            demog_group_def, DemogTblCols, param_names);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%% HIV DIAGNOSIS STATUS TRANSITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%   
        
        %%%%%%%%%%%%%%%%%%%%%%%%%% PREP STATUS TRANSITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%% ADHERENCE STATUS TRANSITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%% ACQUIRING INFECTION %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Derive probability of infection
        
        % Find total people alive (N)
        alive = find_indices(state_matrix, state_mat_demog_group_idx, StateMatCols.alive, '=', 1);
        num_alive = size(alive, 2);
        
        % Find all alive, healthy people (S)
        healthy = find_indices(state_matrix, alive, StateMatCols.hiv, '=', 0);
        num_healthy = size(healthy, 2);
        
        % Find all alive, infected people (I)
        infected = find_indices(state_matrix, alive, StateMatCols.hiv, '>', 0);
        num_infected = size(infected, 2);
        
        
        
        beta = 1
        prob_infection = beta * (num_infected/num_alive) * num_healthy;
        
        

    end   
    
    % Record the state matrix at each point in time
end

