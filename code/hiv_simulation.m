% Run the simulation parameters
hiv_simulation_parameters;

% Read in the initial population from a csv
% Assumption: the csv file has all the people pre-populated with
% characteristics
[state_matrix, StateMatCols] = read_table(init_pop_file);

% Read in all possible demographic groups
[demog_table, DemogTblCols] = create_demog_groups(demog_var_def_file);

% Read in demographic groups specified in the mixing matrix
[mixing_table, MixingTblCols] = create_demog_groups(mixing_mat_def_file);

% Read in the mixing matrix
mixing_matrix = importdata(mixing_mat_file);


% Create array for case when everyone is eligible for subsetting
all_eligible = 1:size(state_matrix, 1);


% For each time period
for t = 1:T
    
    % For each demographic group
    for demog_group_def = demog_table.'
        
        % Get row indices for people in the state matrix in this
        % demographic group
        state_mat_demog_group_idx = find_demog_rows(state_matrix, StateMatCols, demog_group_def, DemogTblCols, DemogTblCols);

        
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
        state_matrix = transition(art_transition_definition_path, art_transition_probabilities_path, ...
                            state_matrix, StateMatCols, state_mat_demog_group_idx, ...
                            demog_group_def, DemogTblCols, param_names);
                        
        %%%%%%%%%%%%%%%%%%%%%%%%%% PREP STATUS TRANSITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%
        state_matrix = transition(art_transition_definition_path, art_transition_probabilities_path, ...
                            state_matrix, StateMatCols, state_mat_demog_group_idx, ...
                            demog_group_def, DemogTblCols, param_names);
                        
        %%%%%%%%%%%%%%%%%%%%%%%%%% ADHERENCE STATUS TRANSITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%
        state_matrix = transition(art_transition_definition_path, art_transition_probabilities_path, ...
                            state_matrix, StateMatCols, state_mat_demog_group_idx, ...
                            demog_group_def, DemogTblCols, param_names);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%% ACQUIRING INFECTION %%%%%%%%%%%%%%%%%%%%%%%%%%
        state_matrix = infection(demog_group_def, DemogTblCols, mixing_matrix, mixing_table, MixingTblCols, state_matrix, StateMatCols);

             

    end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% AGING %%%%%%%%%%%%%%%%%%%%%%%%%%

    % Find everyone who is alive
    alive = find_indices(state_matrix, all_eligible, StateMatCols.alive, '=', 1);

    % Add a unit of time to their ages
    state_matrix(:,StateMatCols.age) = state_matrix(:,StateMatCols.age) + 1;
        
    % Record the state matrix at each point in time
end

