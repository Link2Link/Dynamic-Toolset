function idx_list = get_matched_idx_list_in_cell(cell_data,query_data)
%
% Find the index whose corresponding item matches the 'query_data'
%
idx_list = [];
n = length(cell_data);
for i = 1:n
    if isequal(cell_data{i}, query_data)
        idx_list = [idx_list,i];
    end
end
