function [demog_table, DemogTblCols] = create_demog_groups(demog_var_def_file)
    
    % Read in csv file that defines demographic variable categories
    [demog_var_def, DemogTblCols] = read_table(demog_var_def_file);
        
    % Create possible age_min, age_max pairs
    age_buckets = [demog_var_def(:,DemogTblCols.age_min), demog_var_def(:,DemogTblCols.age_max)];
    
    % Template of variables we've already created combinations for
    template = age_buckets;

    % Placeholder matrix where we can build on top of the template 
    placeholder = [];
    
    % For every variable that isn't age
    for variable = fieldnames(DemogTblCols)'
        variable_string = string(variable);

        if ~strcmp(variable_string, 'age_min') && ~strcmp(variable_string, 'age_max')

            column = demog_var_def(:,DemogTblCols.(variable_string));
            column = column(~isnan(column));
            
            % Number of times to repeat the category
            num_repeats = size(template, 1);
            
            % For every row in each column
            for categ = column'
                to_append = [template repmat(categ, num_repeats, 1)];
                placeholder = [placeholder; to_append]; %#ok<AGROW>
            end 
            
            template = placeholder;
            placeholder = [];
        end
    end
    
    demog_table = template;

end
