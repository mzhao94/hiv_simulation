% Run the simulation parameters
hiv_simulation_parameters;

% Read in the initial population from a csv
% Assumption: the csv file has all the people pre-populated with
% characteristics
% Note: both the alive and in_simulation variables for each person will
% be set to 0 because nobody has been "born" yet and nobody has entered the
% simulation yet

[state_matrix, feature_names] = xlsread( input_file );

% Assign column indices to column names
col_idx = struct() ;
 for i = 1 : length(feature_names)
    col = char(feature_names(1,i));
    col_idx.(col) = i;
 end

for t = 1:T
    
    %%%%%%%%%%%%%%% BIRTHS %%%%%%%%%%%%%%%%%
    
    % Find the people who are eligible to be born (not alive, and not in
    % the simulation)
    eligible_to_be_born = find(state_matrix(:,col_idx.alive) == 0 ...
                                & state_matrix(:,col_idx.in_simulation) == 0);
    
    % TEMP ALGORITHM: Select a certain percentage of people based on a birth_rate
    
    % Don't need to check in_simulation status because if alive == 1, then
    % in_simulation must equal one (the reverse is not necessarily true if
    % a person entered the simulation but died)
    total_alive = length(find(state_matrix(:,col_idx.alive) == 1));
    num_to_be_born = round(birth_rate * total_alive);
    to_be_born = datasample(eligible_to_be_born, num_to_be_born);
    
    % Change these people's state to alive and in_simulation
    state_matrix(to_be_born, col_idx.alive) = 1;
    state_matrix(to_be_born, col_idx.in_simulation) = 1;
    
    
    %%%%%%%%%%%%%%% DEATHS %%%%%%%%%%%%%%%
    
    % Find the people who are eligible to die (alive, and in the
    % simulation)
    eligible_to_die = find(state_matrix(:, col_idx.alive) == 1 ...
                            & state_matrix(:, col_idx.in_simulation) == 1)
    
    % TEMP ALGORITHM: Select a certain percentage of people based on
    % death_rate
    total_alive = length(eligible_to_die);
    num_to_die = round(death_rate * total_alive);
    to_die = datasample(eligible_to_die, num_to_die);
    
    % Change these people's state to not alive, but keep in_simulation to
    % indicate this person has already been in the simulation
    state_matrix(to_die, col_idx.alive) = 0;
    
    %%%%%%%%%%%%%%% IDU TREATMENT ADHERENCE/NON_ADHERENCE %%%%%%%%%%%%%%%
    
    %%% People who start MMT
    
    %%%%%%%%%%%%%%% HIV PREVENTION ADHERENCE/NON-ADHERENCE %%%%%%%%%%%%%%%

    %%% People who start PrEP
    
    % TEMP ALGORITHM: Select a certain percentage of people based on PrEP
    % enrollment rates
    
    % Find people who are alive, are HIV-, and aren't on PrEP
    eligible_for_PrEP = find(state_matrix(:, col_idx.alive) == 1 ...
                            & state_matrix(:, col_idx.HIV) == 0 ...
                            & state_matrix(:, col_idx.PrEP) == 0);
    
    num_PrEP_eligible = length(eligible_for_PrEP);
    num_to_get_PrEP = round(PrEP_rate * num_PrEP_eligible);
    to_get_PrEP = datasample(eligible_for_PrEP, num_to_get_PrEP);
    
    % Change these people's PrEP states
    state_matrix(to_get_PrEP, col_idx.PrEP) = 1;
    
    
    %%% People who were on PrEP but are no longer on it
    
    
    % TEMP ALGORITHM: Select a certain percentage of people based on PrEP
    % enrollment rates
    
    % Find people who are alive, and are on PrEP
    on_PrEP = find(state_matrix(:, col_idx.alive) == 1 ...
                            & state_matrix(:, col_idx.PrEP) == 1);
    
    num_on_PrEP = length(on_PrEP);
    num_to_stop_PrEP = round(PrEP_stop_rate * num_on_PrEP);
    to_stop_PrEP = datasample(on_PrEP, num_to_stop_PrEP);
    
    % Change these people's PrEP states
    state_matrix(to_get_PrEP, col_idx.PrEP) = 0;
    
    
    %%%%%%%%%%%%%%% HEALTH TRANSITIONS %%%%%%%%%%%%%%%
    
    % Assumption: everyone starts off with no treatment, and no HIV
    % TEMP ALGORITHM:
    

    % Determine how many of those people in infected health states get
    % treatment
    
    % Determine how many healthy people acquire some sort of infection
    
    %%%%%%%%%%%%%%% TREATMENT %%%%%%%%%%%%%%%
    
    % Record the state matrix at each point in time
end

