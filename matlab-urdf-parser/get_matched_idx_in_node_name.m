function idx = get_matched_idx_in_node_name(node,name)
%
% Get index of the node matching name
%

idx = 0;
for i = 1:length(node)
    if strcmp(node(i).name,name)
        idx = i;
        return;
    end
end
