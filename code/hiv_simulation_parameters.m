% Number of months to run the simulation
T = 1

% Initial population state matrix input file
input_file = 'input/init_pop.xlsx'

% Birth rate for the population (0.0135 based on current LA county stats)
% Currently arbitrarily high so we can test with small population numbers
birth_rate = 0.5

% Death rate for population
death_rate = 0.1

% The percentage of people without HIV who get PrEP
PrEP_rate = 0.2

% The percentage of people on PrEP that stop
PrEP_stop_rate = 0.1