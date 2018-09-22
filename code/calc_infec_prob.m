% This function will take in a state matrix and calculate the infection
% probabilities for each demographic group using the formula 
% P(infection) = (I/N) * S

function infection_probability = calc_infec_prob(state_matrix, StateMatCols, demog_group_def, demog_col_struct)

    % Find the people in the demographic group
    demog_mat_indices = find_demog_rows(state_matrix, StateMatCols, demog_group_def, demog_col_struct);
    
    % find total number of people in the demographic group (N)
    N = size(demog_mat_indices, 2);

    % calculate the healthy people in the demographic group (S)
    healthy = find_indices(state_matrix, demog_mat_indices, StateMatCols.hiv, '=', 0);
    S = size(healthy, 2);

    % find all the transmitting people in the demogrphic group (I) -
    % defined by their hiv status of at least 1, and prep and art status of 0
    has_hiv = find_indices(state_matrix, demog_mat_indices, StateMatCols.hiv, '>', 0);
    has_hiv_no_art = find_indices(state_matrix, has_hiv, StateMatCols.art, '=', 0);
    has_hiv_no_art_no_prep = find_indices(state_matrix, has_hiv_no_art, StateMatCols.prep, '=', 0);
    I = size(has_hiv_no_art_no_prep, 2);
    
    infection_probability = (I/N) * S;

end