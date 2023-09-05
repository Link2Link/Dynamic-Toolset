ccc
%% Animiate random trajectory
ccc

% Save video
MAKE_VIDEO = 0;

% Load the model
model_name = 'sawyer'; % franka / sawyer
urdf_path = sprintf('model/%s/%s_urdf.xml',model_name,model_name);
base_name = 'world_fixed'; VERBOSE = 0;
model = parse_urdf(model_name,urdf_path,base_name,VERBOSE);

% Initialize joints
names2control = model.joint_names_revolute;
positions2control_rad = zeros(1,length(names2control));
init_units = p2u_urdf(model,names2control,positions2control_rad);
model = set_positions_urdf(model,names2control,positions2control_rad,VERBOSE);

% Get GRP
t_min = 0; t_max = 1; n_anchor = 2; n_test = 100;
t_anchor = linspace(t_min,t_max,n_anchor)';
x_anchor = [init_units';init_units'];
t_test = linspace(t_min,t_max,n_test)';
l_anchor = 1.0*ones(n_anchor,1); l_test = 1.0*ones(n_test,1);
kfun_str = 'kernel_levrq';
hyp_mu = [0.1,0.1,1.0]; % [beta/gamma/alpha]
sig2w_mu = 1e-8;
hyp_var = [0.2,0.1,1.0]; % [beta/gamma/alpha]
sig2w_var = 1e-8;
grp = get_grp(t_anchor,x_anchor,t_test,l_anchor,l_test,...
    kfun_str,hyp_mu,sig2w_mu,hyp_var,sig2w_var);
% Plot GRP
opt = struct('n_sample',1,'seed',0,'names','','title_str','');
plot_grp(grp,opt);

%%
% Sample from GRP
n_sample = 1;
sampled_path_list = sample_grp(grp,n_sample);
sampled_path = sampled_path_list{1};

% Plot random positions
opt = struct('fig_idx',2,'subfig_idx',1,...
    'view_info',[116,17],'fig_size',[0.4,0.2,0.55,0.7],...
    'ball_r',0.02,'coord_r',0.05,...
    'title_str',model_name,'title_fs',25,...
    'axis_info',[-1.5,+1.5,-1.5,+1.5,0.0,1.8],...
    'SHOW_TEXT',0,'PLOT_LINKS',0,'PLOT_JOINTS',0,'PLOT_COORDINATES',0, ...
    'PLOT_ROTATIONAL_AXES',0,'PLOT_MESH',1,'PLOT_CUBE',1);
for tick = 1:n_test % for all test time steps
    init_units = sampled_path(tick,:);
    model = set_positions_urdf(model,names2control,u2p_urdf(model,names2control,init_units),VERBOSE);
    opt.title_str = sprintf('[%d/%d] %s',tick,n_test,model_name);
    fig = plot_model_urdf(model,opt);
    drawnow;
    if MAKE_VIDEO
        png_name = sprintf('temp/pic_%03d.png',tick);
        save_fig2png(fig,png_name,VERBOSE);
    end
end % for i = 1:n_test % for all test time steps

% Save video 
if MAKE_VIDEO 
    curr_time = datestr(now,'yy-mm-dd_HH-MM-SS');
    vid_name = sprintf('vid/vid_grp_%s_%s.mp4',...
        model_name,curr_time);
    frm_rate = 20;
    video = VideoWriter( vid_name, 'MPEG-4' );
    video.FrameRate = ( frm_rate );
    open( video );
    fprintf('Start saving a video..');
    for tick = 1:n_test % For all time (sec)
        loadPicName = sprintf('temp/pic_%03d.png',tick);
        img = imread( loadPicName ); img = im2double(img);
        writeVideo( video, img );
    end
    close( video );
    fprintf('[%s] saved.\n',vid_name);
end

%%
