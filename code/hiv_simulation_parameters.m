% Number of months to run the simulation
T = 1;

% Initial population state matrix input file
init_pop_file = 'input/init_pop.csv';

% File path for table that defines all demographic variable categories
demog_var_def_file = 'input/demog_var_def.csv';

% Paths defined for birth transitions
birth_transition_definition_path = 'input/birth_transitions/state_definitions/';
birth_transition_probabilities_path = 'input/birth_transitions/state_probabilities/';

% Paths defined for death transitions
death_transition_definition_path = 'input/death_transitions/state_definitions/';
death_transition_probabilities_path = 'input/death_transitions/state_probabilities/';

% Paths defined for hiv status transitions
hiv_transition_definition_path = 'input/hiv_transitions/state_definitions/';
hiv_transition_probabilities_path = 'input/hiv_transitions/state_probabilities/';

% Paths defined for art status transitions
art_transition_definition_path = 'input/art_transitions/state_definitions/';
art_transition_probabilities_path = 'input/art_transitions/state_probabilities/';


% Define distribution parameter names
param_names = {'alpha', 'beta'};



%%% Define categorical variable values as explicit names for easy
%%% interpretation

% Race
num_race_categ = 5;
white = 0;
black = 1;
latinX = 2;
asian_pac = 3;
other = 4;

% HIV
num_hiv_categ = 4;
no_hiv = 0;
early_hiv = 1;
late_hiv = 2;
aids = 3;

% PrEP
num_prep_categ = 3;
no_prep = 0;
prep_adh = 1;
prep_non_adh = 2;

% ART
num_art_categ = 3;
no_art = 0;
art_adh = 1;
art_non_adh = 2;






