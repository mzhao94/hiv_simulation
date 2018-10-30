hiv_simulation;

results = prevalence_by_time(state_matrices_path, StateMatCols, 'hiv', [1, 0], '','')

plot(1:T, results);