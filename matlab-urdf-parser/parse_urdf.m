function model = parse_urdf(model_name,urdf_path,base_name,VERBOSE)
%
% Parse URDF file
%

% Load the URDF file
s = urdf2struct(urdf_path);

% Parse all the information
robot = s.robot;
model.name = model_name;
model.robot = robot;
model.n_joint = numel(robot.joint);
model.n_link = numel(robot.link);

model.joint_names = cell(model.n_joint,1);
model.joint_names_revolute = ''; temp_idx = 0;
model.parent_links = cell(model.n_joint,1);
model.child_links = cell(model.n_joint,1);

model.joint_rpys = cell(model.n_joint,1);
model.joint_xyzs = cell(model.n_joint,1);
model.joint_axes = cell(model.n_joint,1);
model.joint_limits = cell(model.n_joint,1);
model.joint_types = cell(model.n_joint,1);

for i = 1:model.n_joint % for all joints
    robot_joint_i = struct2num(robot.joint{i});
    % parse information
    model.urdf(i).joint_name = robot_joint_i.Attributes.name; % joint name
    model.urdf(i).type = robot_joint_i.Attributes.type; % movable or not
    if isequal(model.urdf(i).type,'revolute')
        temp_idx = temp_idx + 1;
        model.joint_names_revolute{temp_idx} = model.urdf(i).joint_name;
    end
    model.urdf(i).parent_link = robot_joint_i.parent.Attributes.link; % parent link
    model.urdf(i).child_link = robot_joint_i.child.Attributes.link; % child link
    % joint coordinate
    model.urdf(i).joint_rpy = robot_joint_i.origin.Attributes.rpy;
    model.urdf(i).joint_xyz = robot_joint_i.origin.Attributes.xyz;
    % some joints do not have field: axis
    try
        model.urdf(i).joint_axis = robot_joint_i.axis.Attributes.xyz;
    catch
        model.urdf(i).joint_axis = [0,0,0];
    end
    % some joints do not have field: limit
    try
        joint_limit = robot_joint_i.limit.Attributes;
        range = [joint_limit.lower,joint_limit.upper];
        model.urdf(i).joint_limit = range;
    catch
        model.urdf(i).joint_limit = [0,0];
    end
    % Append joint and link names for further search
    model.joint_names{i} = model.urdf(i).joint_name;
    model.parent_links{i} = model.urdf(i).parent_link;
    model.child_links{i} = model.urdf(i).child_link;
    model.joint_rpys{i} = model.urdf(i).joint_rpy;
    model.joint_xyzs{i} = model.urdf(i).joint_xyz;
    model.joint_axes{i} = model.urdf(i).joint_axis;
    model.joint_limits{i} = model.urdf(i).joint_limit;
    model.joint_types{i} = model.urdf(i).type;
end % for i = 1:model.n_joint % for all joints

% Parse link information
model.link_rpys = cell(model.n_link,1);
model.link_xyzs = cell(model.n_link,1);
model.link_filenames = cell(model.n_link,1);
model.link_scales = cell(model.n_link,1);
model.link_names = cell(model.n_link,1);

model.link_fvs = cell(model.n_link,1);
model.link_scaled_vertices = cell(model.n_link,1);
model.link_vmins = cell(model.n_link,1);
model.link_vmaxs = cell(model.n_link,1);
model.link_vlens = cell(model.n_link,1);

model.link_ps = cell(model.n_link,1);
model.link_Rs = cell(model.n_link,1);

for i = 1:model.n_link % for all links
    robot_link_i = robot.link{i};
    model.link_names{i} = robot_link_i.Attributes.name;
    
    if isfield(robot_link_i,'visual')
        if isfield(robot_link_i.visual,'origin')
            model.link_rpys{i} = str2num(robot_link_i.visual.origin.Attributes.rpy);
            model.link_xyzs{i} = str2num(robot_link_i.visual.origin.Attributes.xyz);
        else
            model.link_rpys{i} = [0,0,0];
            model.link_xyzs{i} = [0,0,0];
        end
        try
            model.link_scales{i} = str2num(robot_link_i.visual.geometry.mesh.Attributes.scale);
        catch
            model.link_scales{i} = [1,1,1];
        end
        try
            [~,name,ext] = fileparts(robot_link_i.visual.geometry.mesh.Attributes.filename);
            if ~isequal(ext,'.stl')
                ext = '.stl';
            end
            model.link_filenames{i} = [name,ext];
        catch
            model.link_filenames{i} = '';
        end
    else
        model.link_rpys{i} = [0,0,0];
        model.link_xyzs{i} = [0,0,0];
        model.link_filenames{i} = '';
        model.link_scales{i} = [1,1,1];
    end
    
    % Load STL if exists
    if model.link_filenames{i}
        fv = stlread(['model/',model.name,'/visual/',model.link_filenames{i}]); % read stl
        n_vertice = size(fv.vertices,1);
        reduce_patch_ratio = min(1,floor(100000/n_vertice));
        reduce_patch_ratio = max(0.1,reduce_patch_ratio);
        fv = reducepatch(fv,reduce_patch_ratio); % reduce the resolution
        fv.vertices = model.link_scales{i}.*fv.vertices; % scale the mesh
        v_min = min(fv.vertices);
        v_max = max(fv.vertices);
        v_len = v_max-v_min;
        % append
        model.link_fvs{i} = fv;
        model.link_scaled_vertices{i} = fv.vertices;
        model.link_vmins{i} = v_min;
        model.link_vmaxs{i} = v_max;
        model.link_vlens{i} = v_len;
    end
end

% print-out information (joint, links, offsets ...)
if VERBOSE
    for i = 1:model.n_joint % for all joints
        str_temp = ['[%02d/%d] joint:[%12s] (type:%9s) parent:[%18s] '...
            'child:[%18s] rpy_offset:%s xyz:%s axis:%s limit:%s \n'];
        if isequal(model.urdf(i).type,'fixed')
            print_col = 2;
            fprintf(print_col,str_temp , ...
                i,model.n_joint,model.urdf(i).joint_name,...
                model.urdf(i).type,...
                model.urdf(i).parent_link,...
                model.urdf(i).child_link,...
                vec2str(model.urdf(i).joint_rpy,'%+.2f'),...
                vec2str(model.urdf(i).joint_xyz,'%+.2f'),...
                vec2str(model.urdf(i).joint_axis,'%+d'),...
                vec2str(model.urdf(i).joint_limit,'%+.2f') ...
                );
        else
            print_col = 1;
            fprintf(print_col,str_temp , ...
                i,model.n_joint,model.urdf(i).joint_name,...
                model.urdf(i).type,...
                model.urdf(i).parent_link,...
                model.urdf(i).child_link,...
                vec2str(model.urdf(i).joint_rpy,'%+.2f'),...
                vec2str(model.urdf(i).joint_xyz,'%+.2f'),...
                vec2str(model.urdf(i).joint_axis,'%+d'),...
                vec2str(model.urdf(i).joint_limit,'%+.2f') ...
                );
        end
    end
end

% Get unique joint names
model.unique_joint_names = cell(model.n_joint,1);
for i = 1:model.n_joint % for all joints
    model.unique_joint_names{i} = model.urdf(i).joint_name;
end

% Get unique link names
idx = 0;
model.unique_link_names = cell(2,1);
for i = 1:model.n_joint % for all joints
    if isempty(get_matched_idx_list_in_cell(model.unique_link_names,model.urdf(i).parent_link))
        idx = idx + 1;
        model.unique_link_names{idx} = model.urdf(i).parent_link;
    end
    if isempty(get_matched_idx_list_in_cell(model.unique_link_names,model.urdf(i).child_link))
        idx = idx + 1;
        model.unique_link_names{idx} = model.urdf(i).child_link;
    end
end

% Make a tree structure ('name','childs') from mode.urdf
% => 'model.node' will have 'name' and 'childs' as its member
for i = 1:model.n_joint
    model.node(i) = struct('name',model.joint_names{i},'childs',[]); % set joint names
    childs = get_matched_idx_list_in_cell(model.parent_links,model.child_links{i});
    model.node(i).childs = childs;
end

% Update the name of the base node
model.base_name = base_name;

% Add parents to node
base_idx_in_node = get_matched_idx_in_node_name(model.node,model.base_name);
model.node = add_parent2node(model.node,base_idx_in_node);

% Convert to LCRS tree (total number of nodes remains the same)
model.node = convert_multiwayTree_to_lcrsTree(model.node,base_name);

% Initialize kinematic chain
model = init_kinematic_chain_for_urdf(model);

% Add positions in radian
model.positions_rad = zeros(model.n_joint,1); % set positions in radian


