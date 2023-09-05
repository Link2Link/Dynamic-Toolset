function node = convert_multiwayTree_to_lcrsTree(node,base_name)
%
% Convert an arbitrary multiway tree to LCRS tree
%  (Left Child / Right Sibling)
%j

% Add depth to node
base_idx = get_matched_idx_in_node_name(node,base_name);
node = add_depth2node(node,base_idx);

% Initialize LCRS structure
for i = 1:length(node)
    node(i).lcrsChild = 0;
    node(i).lcrsSibling = 0;
end

% Compute depth-index list ('depth2idxList')
maxDepth = 0;
for i = 1:length(node) % compute max depth 
    if node(i).depth > maxDepth, maxDepth=node(i).depth; end
end
depth2idxList = cell(maxDepth,1);
for d = 1:maxDepth
    depth2idxList{d} = [];
end
for i = 1:length(node)
    d = node(i).depth;
    depth2idxList{d} = [depth2idxList{d}, i];
end

% Loop over depth
for d = 1:maxDepth % For each depth
    for i = depth2idxList{d} % For each index in current depth 
        if d == 1
            % Root node
        else
            % Other nodes
            % 1. Find parent's node
            parentIdx = node(i).parent;
            
            % Add current node to parent or sibling
            if node(parentIdx).lcrsChild == 0
                % If parent's child does not exist, add current node to parent
                node(parentIdx).lcrsChild = i;
            else
                % If parent already has child, add current node to parent's child's sibling
                childIdx = node(parentIdx).lcrsChild;
                node = add_sibling(node,childIdx,i); % <= Add! 
            end
        end
    end
end

% Add 'lcrsParent' to nodes
base_idx = get_matched_idx_in_node_name(node,base_name);
node = add_lcrsParent2node(node,base_idx);



% ------------------------------------------------------------------- %
function node = add_sibling(node,childIdx,idx)
%
% Add sibling to the node with recursion 
% 
if node(childIdx).lcrsSibling == 0
    node(childIdx).lcrsSibling = idx;
    return;
else
    % Recursive 
    node = add_sibling(node,node(childIdx).lcrsSibling,idx);
end
% ------------------------------------------------------------------- %

% ------------------------------------------------------------------- %
function node = add_lcrsParent2node(node,idx)
%
% Add LCRS parent to nodes with recursion
%
if nargin == 1
    idx = 1;
end
lcrsChildIdx = node(idx).lcrsChild;
if lcrsChildIdx == 0
    return;
end

% Add parent 
node(lcrsChildIdx).lcrsParent = idx;

% Recursive
node = add_lcrsParent2node(node,lcrsChildIdx);
% ------------------------------------------------------------------- %