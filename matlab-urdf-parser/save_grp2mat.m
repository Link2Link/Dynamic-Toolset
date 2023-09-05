function save_grp2mat(grp,mat_name)
%
% Save GRP to a mat file
%

% Parse 
t_anchor = grp.t_anchor;
x_anchor = grp.x_anchor;
t_test = grp.t_test;
l_anchor = grp.l_anchor;
l_test = grp.l_test;
kfun_str = grp.kfun_str;
hyp_mu = grp.hyp_mu;
sig2w_mu = grp.sig2w_mu;
hyp_var = grp.hyp_var;
sig2w_var = grp.sig2w_var;

save(mat_name,'t_anchor','x_anchor','t_test','l_anchor','l_test',...
    'kfun_str','hyp_mu','sig2w_mu','hyp_var','sig2w_var');
fprintf('[%s] saved.\n',mat_name);
