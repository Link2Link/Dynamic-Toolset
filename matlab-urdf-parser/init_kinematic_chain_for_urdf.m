function model = init_kinematic_chain_for_urdf(model,idx_to)
%
% Initialize the kinematic chain for URDF files
%

if nargin == 1
    idx_to = 1;
end
if idx_to == 0, return; end

% T and R will be updated from FK not here. 
model.node(idx_to).T = zeros(3,1);
model.node(idx_to).R = eye(3,3);

% b, axis, and rpy_offset are paresed from the model
model.node(idx_to).b = model.joint_xyzs{idx_to}';
model.node(idx_to).axis = model.joint_axes{idx_to}; % signed rotational axis 
model.node(idx_to).rpy_rad_offset = model.joint_rpys{idx_to}; % rpy offset 

% rpy will be updated with 'model.positions_rad' combined with 'model.node(idx_to).axis'
% model.node(idx_to).position_rad = 0;
model.node(idx_to).rpy_rad = zeros(3,1); % this will simply be updated. 

% Position range
model.node(idx_to).limit = model.joint_limits{idx_to};

% Recursive
model = init_kinematic_chain_for_urdf(model,model.node(idx_to).lcrsSibling);
model = init_kinematic_chain_for_urdf(model,model.node(idx_to).lcrsChild);
