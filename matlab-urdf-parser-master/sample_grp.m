function sampled_path_list = sample_grp(grp,n_sample)
%
% Sample paths from the GRP distribution. 
%

sampled_path_list = cell(n_sample,1);
for i = 1:n_sample % for each samples
    sampled_path = zeros(grp.n_test,grp.xdim);
    for d_idx = 1:grp.xdim % for each dimension
        % R_mtx = randn(grp.n_test,1);
        % sampled_path(:,d_idx) = chol(grp.K)*R_mtx;
        sampled_path(:,d_idx) = mvnrnd(zeros(grp.n_test,1),grp.K,1)';
    end
    sampled_path_list{i} = sampled_path+grp.mu;
end
