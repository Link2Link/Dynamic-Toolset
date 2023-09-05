function model = forward_kinematics_urdf(model,idx_to)
%
% Forward kinematics the kinematic chain of the model.node %
% -> ('R'/'T'/'roll'/'pitch'/'b')
%
% Here, we use radian for angles.
%

if nargin == 1
    idx_to = get_matched_idx_in_node_name(model.node,model.base_name);
end
if idx_to == 0, return; end

% Do things here
idx_fr = model.node(idx_to).parent;

% Forward Kinematics
if ~isempty(idx_fr)
    % If the parent model.node exists, do FK:
    % Do translation first
    model.node(idx_to).T = model.node(idx_fr).R*model.node(idx_to).b+model.node(idx_fr).T;
    % Compute rpy
    model.node(idx_to).rpy_rad = model.positions_rad(idx_to)*model.joint_axes{idx_to};
    rpy_deg = model.node(idx_to).rpy_rad*180/pi;
    R_ctrl = get_rotMtx_rpy_deg(rpy_deg(1),rpy_deg(2),rpy_deg(3));
    rpy_deg = model.node(idx_to).rpy_rad_offset*180/pi;
    R_offset = get_rotMtx_rpy_deg(rpy_deg(1),rpy_deg(2),rpy_deg(3));
    model.node(idx_to).R = model.node(idx_fr).R*R_offset*R_ctrl; % offset first then control
else
    DO_NOTHING = true;
end

% Recursive
model = forward_kinematics_urdf(model,model.node(idx_to).lcrsSibling);
model = forward_kinematics_urdf(model,model.node(idx_to).lcrsChild);
