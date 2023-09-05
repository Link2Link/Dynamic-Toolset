ccc
%% Parse URDF
ccc

% Load the model
model_name = 'sawyer'; % franka / sawyer
urdf_path = sprintf('model/%s/%s_urdf.xml',model_name,model_name);
base_name = 'world_fixed'; VERBOSE = true;
model = parse_urdf(model_name,urdf_path,base_name,VERBOSE);

% Plot parsed URDF tree structure
opt = struct('fig_idx',1,'fig_size',[0.1,0.4,0.8,0.3],...
    'title_str','Parsed from URDF');
fig_node = plot_node(model.node,opt);

% Initialize joints
names2control = model.joint_names_revolute;
positions2control_rad = zeros(1,length(names2control));
model = set_positions_urdf(model,names2control,positions2control_rad,VERBOSE);

% Plot the kinematic chain with mesh
opt = struct('fig_idx',2,'subfig_idx',1,...
    'view_info',[120,24],'fig_size',[0.4,0.2,0.55,0.7],...
    'ball_r',0.02,'coord_r',0.05,...
    'title_str',model_name,'title_fs',25,...
    'axis_info',[-1.5,+1.5,-1.5,+1.5,0.0,1.6],...
    'SHOW_TEXT',1,'PLOT_LINKS',1,'PLOT_JOINTS',0,'PLOT_COORDINATES',1, ...
    'PLOT_ROTATIONAL_AXES',1,'PLOT_MESH',1,'PLOT_CUBE',1);
fig_model = plot_model_urdf(model,opt);
drawnow;

% Plot the simple kinematic chain
opt = struct('fig_idx',3,'subfig_idx',1,...
    'view_info',[120,24],'fig_size',[0.4,0.1,0.55,0.7],...
    'ball_r',0.1,'coord_r',0.1,...
    'title_str',model_name,'title_fs',25,...
    'axis_info',[-1.5,+1.5,-1.5,+1.5,0.0,1.6],...
    'SHOW_TEXT',1,'PLOT_LINKS',1,'PLOT_JOINTS',1,'PLOT_COORDINATES',1, ...
    'PLOT_ROTATIONAL_AXES',1,'PLOT_MESH',0,'PLOT_CUBE',0);
fig_simple = plot_model_urdf(model,opt);
drawnow;


% Save figure
png_name = sprintf('fig/fig_node_%s.png',model_name);
save_fig2png(fig_node,png_name,VERBOSE);
png_name = sprintf('fig/fig_model_%s.png',model_name);
save_fig2png(fig_model,png_name,VERBOSE);
png_name = sprintf('fig/fig_simple_%s.png',model_name);
save_fig2png(fig_simple,png_name,VERBOSE);

%% Plot random positions
max_tick = 10;
for tick = 1:max_tick
    units = rand(length(names2control),1); % random joint units (0~1)
    model = set_positions_urdf(model,names2control,u2p_urdf(model,names2control,units),VERBOSE);
    opt.title_str = sprintf('[%d/%d] %s',tick,max_tick,model_name);
    plot_model_urdf(model,opt);
    drawnow;
end

%% 
