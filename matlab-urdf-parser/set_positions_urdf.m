function model = set_positions_urdf(model,names,positions_rad,VERBOSE)
%
% Set positions to the model from URDF 
%

for i = 1:length(names)
    % find the index
    name = names{i};
    idx = get_matched_idx_in_node_name(model.node,name);
    if idx == 0
        fprintf(2,'[set_positions_urdf] cannot find [%s].\n',name);
    else
        if ~isequal(model.joint_types{idx},'revolute')
            fprintf(2,'[set_positions_urdf] [%s] cannot move.\n',name);
            return;
        end
        model.positions_rad(idx) = positions_rad(i); % update positions 
    end
end

% Run FK
model = trim_positions_urdf(model,VERBOSE);
model = forward_kinematics_urdf(model);

% Update mesh
model = update_mesh_urdf(model);
