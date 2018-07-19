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
    col = char(feature_names(1,i))
    col_idx.(col) = i;
 end

for t = 1:T
    
    %%%%%%%%%%%%%%% BIRTHS %%%%%%%%%%%%%%%%%
    
    % Find the people who are eligible to be born (not alive, and not in
    % the simulation)
    eligible_to_be_born = find(state_matrix(:,col_idx.alive) == 0 & state_matrix(:,col_idx.in_simulation) == 0);
    
    % TEMP ALGORITHM: Select a certain percentage of people based on a birth_rate
    to_be_born = datasample(data,k)
    
    % Determine how many of those people that are born die
    
    
    % Determine how many of those people transition to/from each health state 


    % Determine how many of those people in infected health states get
    % treatment
    
    % Determine how many healthy people acquire some sort of infection
    
    % Record the state matrix at each point in time
end

