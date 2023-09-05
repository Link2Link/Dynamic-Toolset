ccc
%% Gaussian random path
ccc

% Set data for GRP
n_anchor = 2; dim = 3; n_test = 100;
t_min = 0; t_max = 1;
t_anchor = linspace(t_min,t_max,n_anchor)';
x_anchor = rand(n_anchor,dim);
t_test = linspace(t_min,t_max,n_test)';
l_anchor = 1.0*ones(n_anchor,1);
l_test = 1.0*ones(n_test,1);
kfun_str = 'kernel_levrq';
hyp_mu = [1,0.1,1.0]; % [beta/gamma/alpha]
sig2w_mu = 1e-8;
hyp_var = [1,0.1,1.0]; % [beta/gamma/alpha]
sig2w_var = 1e-8;
% Get GRP
grp = get_grp(t_anchor,x_anchor,t_test,l_anchor,l_test,...
    kfun_str,hyp_mu,sig2w_mu,hyp_var,sig2w_var);

% Plot GRP
opt = struct('n_sample',10,'seed',0,'names','','title_str',''); 
plot_grp(grp,opt);

% Save GRP 
mat_name = 'data/grp_from_matlab.mat';
save_grp2mat(grp,mat_name);

%% Load GRP
ccc;
mat_name = 'data/grp_from_matlab.mat';
grp_loaded = load_grp_from_mat(mat_name);

% Plot loaded GRP 
opt = struct('n_sample',10,'seed',0,'names','','title_str',''); 
plot_grp(grp_loaded,opt);

%%
