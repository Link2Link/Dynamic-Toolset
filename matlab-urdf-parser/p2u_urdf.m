function units = p2u_urdf(model,names,positions)
%
% actual positions => unit positions varying from 0 to 1 
%

n_position = length(names);
units = zeros(n_position,1);

for p_idx = 1:n_position
    curr_name = names{p_idx};
    curr_position = positions(p_idx);
    idx = get_matched_idx_in_node_name(model.node,curr_name);
    joint_limit = model.joint_limits{idx};
    units(p_idx) = (curr_position-joint_limit(1)) / ...
        (joint_limit(2)-joint_limit(1));
end
