function model = update_mesh_urdf(model)
%
% Update mesh information
%
for i = 1:model.n_joint % for all joint
    % get T and R of the current joint
    joint_p = model.node(i).T; joint_R = model.node(i).R;
    % get child link of the current joint
    curr_child_link_name = model.child_links{i};
    curr_child_link_idx = ...
        get_matched_idx_list_in_cell(model.link_names,curr_child_link_name);
    curr_link_rpy = model.link_rpys{curr_child_link_idx};
    curr_link_xyz = model.link_xyzs{curr_child_link_idx}';
    curr_link_filename = model.link_filenames{curr_child_link_idx};
    % if stl exists, load and plot
    if curr_link_filename
        % Get item
        fv = model.link_fvs{curr_child_link_idx};
        scaled_vertices = model.link_scaled_vertices{curr_child_link_idx};
        % Move vertices
        rpy_deg = curr_link_rpy*180/pi; % this is the offset
        p_offset = curr_link_xyz;
        R_offset = get_rotMtx_rpy_deg(rpy_deg(1),rpy_deg(2),rpy_deg(3));
        T_joint = [joint_R, joint_p; 0,0,0,1];
        T_offset = [R_offset, p_offset; 0,0,0,1];
        T_new = T_offset*T_joint; %
        % Get p and R of updated T
        p_link = T_new(1:3,4);
        R_link = T_new(1:3,1:3);
        % Update the mesh 
        V = scaled_vertices*(R_link)' ...
            + ones(size(scaled_vertices,1),1)*(p_link)';
        fv.vertices = V;
        % Append
        model.link_ps{curr_child_link_idx} = p_link;
        model.link_Rs{curr_child_link_idx} = R_link;
        model.link_fvs{curr_child_link_idx} = fv;
    end
end % for i = 1:model.n_joint % for all joint
