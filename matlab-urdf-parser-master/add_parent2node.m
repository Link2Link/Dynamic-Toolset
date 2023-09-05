function node = add_parent2node(node,idx)
%
% Add parent to each node
%
if nargin == 1
    idx = 1;
end
for childIdx = node(idx).childs
    % Add parent to current child
    node(childIdx).parent = idx;
    node = add_parent2node(node,childIdx);
end
