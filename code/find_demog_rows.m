% Get the row corresponding to the demographic group
function demog_mat_indices = find_demog_rows(target_table, target_table_col_struct, demog_group_def, demog_col_struct)
    demog_mat_indices = 1:size(target_table,1);
    for field = fieldnames(demog_col_struct).'
        field_string = string(field);
        target_table_field_col = target_table_col_struct.(field_string);
        demog_field_value = demog_group_def(demog_col_struct.(field_string));
        demog_mat_indices = find_indices(target_table, demog_mat_indices, target_table_field_col, '=', demog_field_value);
    end
end