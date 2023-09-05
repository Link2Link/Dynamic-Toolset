function grp = get_grp(t_anchor,x_anchor,t_test,l_anchor,l_test,...
    kfun_str,hyp_mu,sig2w_mu,hyp_var,sig2w_var)
%
% Get Gaussian Random Path
%

n_anchor = size(t_anchor,1);
n_test = size(t_test,1);
xdim = size(x_anchor,2);
kfun = str2func(kfun_str);

% Make GRP mean zero
x_anchor_mu = mean(x_anchor);
x_anchor_mz = x_anchor - x_anchor_mu;

% Define GRP-mu kernel matrices
l_anchor_mu = ones(size(l_anchor));
K_ta_mu = kfun(t_test,t_anchor,l_test,l_anchor_mu,hyp_mu);
K_aa_mu = kfun(t_anchor,t_anchor,l_anchor_mu,l_anchor_mu,hyp_mu) ...
    + sig2w_mu*eye(n_anchor,n_anchor);
grp_mu = K_ta_mu/K_aa_mu*x_anchor_mz + x_anchor_mu;

% Define GRP-var kernel matrices
K_ta_var = kfun(t_test,t_anchor,l_test,l_anchor,hyp_var);
K_aa_var = kfun(t_anchor,t_anchor,l_anchor,l_anchor,hyp_var) ...
    + sig2w_var*eye(n_anchor,n_anchor);
K_tt_var = kfun(t_test,t_test,l_test,l_test,hyp_var) + 1e-12*eye(n_test,n_test);
grp_var_mtx = K_tt_var - K_ta_var/K_aa_var*K_ta_var' + 1e-12*eye(n_test,n_test);
grp_var_mtx = 0.5*(grp_var_mtx+grp_var_mtx')/2;
grp_var = diag(grp_var_mtx);
grp_std = sqrt(max(0,grp_var));

% Save grp struct 
grp = struct('mu',grp_mu,'var',grp_var,'std',grp_std,...
    't_anchor',t_anchor,'x_anchor',x_anchor,'t_test',t_test,...
    'l_anchor',l_anchor,'l_test',l_test,...
    'n_anchor',n_anchor,'n_test',n_test,...
    'xdim',xdim,...
    'K',grp_var_mtx,'x_anchor_mu',x_anchor_mu,...
    'kfun_str',kfun_str,'kfun',kfun,...
    'hyp_mu',hyp_mu,'sig2w_mu',sig2w_mu,...
    'hyp_var',hyp_var,'sig2w_var',sig2w_var);

%%