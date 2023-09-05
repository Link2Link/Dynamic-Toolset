function positions = u2p_urdf(model,names,units)
%
% unit positions varying from 0 to 1 => actual positions
%

n_position = length(names);
positions = zeros(n_position,1);

for p_idx = 1:n_position
    curr_name = names{p_idx};
    curr_unit = units(p_idx);
    idx = get_matched_idx_in_node_name(model.node,curr_name);
    joint_limit = model.joint_limits{idx};
    positions(p_idx) = joint_limit(1)+...
        (joint_limit(2)-joint_limit(1))*curr_unit;
end
