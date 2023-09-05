function grp = load_grp_from_mat(mat_name)
%
% Load GRP from a mat file
%

l = load(mat_name);
grp = get_grp(l.t_anchor,l.x_anchor,l.t_test,l.l_anchor,l.l_test,...
    l.kfun_str,l.hyp_mu,l.sig2w_mu,l.hyp_var,l.sig2w_var);
