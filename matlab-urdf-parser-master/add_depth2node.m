function node = add_depth2node(node,idx)
%
% Add depth of each node
%

parent_idx = node(idx).parent;
if isempty(parent_idx)
    node(idx).depth = 1;
else
    node(idx).depth = node(parent_idx).depth + 1;
end


% --------------------------------------------------------- %
% Recursion
for child_idx = node(idx).childs
    node = add_depth2node(node,child_idx);
end
