function K = kernel_levrq(X1, X2, L1, L2, hyp)
%
% Rational Quadratic Kernel
%
n1 = size(X1, 1);
n2 = size(X2, 1);
d1 = size(X1, 2);
d2 = size(X2, 2);
if d1 ~= d2, fpinrtf(2, 'Data dimension missmatch! \n'); end

beta  = hyp(1);  % gain
gamma = hyp(2:end-1); % length parameters 
alpha = hyp(end); % RQ alpha

% Make each leverave vector a column vector
L1 = reshape(L1,[],1);
L2 = reshape(L2,[],1);

% Compute Kernel Matrix
x_dists = pdist2(X1./gamma, X2./gamma, 'euclidean').^2;
l_dists = pdist2(L1, L2, 'cityblock');

K = beta*(1+1/2/alpha*x_dists).^(-alpha) ...
    .*cos(pi/2*l_dists) ;

% Limit condition number
if n1 == n2 && n1 > 1 && 0
    sig = 1E-12*beta;
    K = K + sig*eye(size(K));
end
