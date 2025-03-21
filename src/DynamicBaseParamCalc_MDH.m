function base = DynamicBaseParamCalc_MDH(rbt_df, dyn)
    tic;
    format long;
    disp("Calculating base parameter...");
    [r, base_params ,E,beta] = find_dyn_param_deps(rbt_df, dyn); % base on QR
    disp('calc_base_param_: time ' + string(toc) + 'sec');

    n = length(rbt_df.link_params);
    s = n-r;
    % KG = [eye(r), Kd; zeros(s, r), eye(s)];
    % eye_param_num = eye(n);
    % inertia2beta = KG * eye_param_num(:, P)';

    prime_params = rbt_df.link_params';
    % beta_param = inertia2beta * prime_param;
    % base_param = beta_param(1:r);


    H_b = dyn.H * E(:,1:r);

    base.base_num = r;
    base.prime_num = n;
    base.redundant_num = s;
    base.prime_params = prime_params;
    % base.beta_param = beta_param;
    base.base_params = base_params;
    % base.prime2beta = inertia2beta;
    base.H_b = collect(expand(simplify(H_b)));
    base.H = dyn.H;
    
    
    
end


function [r, base_params, E,beta] = find_dyn_param_deps(rbt_df, dyn)
    dof = rbt_df.dof;
    param_num = length(rbt_df.link_params);
    sample_num = param_num*2;
    W = zeros(dof * sample_num, param_num);
    
    funcHandle = matlabFunction(vpa(dyn.H),'Vars', {rbt_df.q_var', rbt_df.dq_var', rbt_df.ddq_var'});
    
    for i = 1:sample_num
        q = rand([dof, 1])*2*pi - pi;
        dq = rand([dof, 1])*2*pi - pi;
        ddq = rand([dof, 1])*2*pi - pi;
        
        W(((i-1)*dof+1):(i*dof), :) = funcHandle(q,dq,ddq);
    end

    %% QR decomposition: W*E = Q*R
    %   R is an upper triangular matrix
    %   Q is a unitary matrix
    %   E is a permutation matrix
    [Q, R, E] = qr(W); 

    % The rank of W is b, which is the number of base parameters
    r = rank(W);
    
    % R = [R1 R2; 
    %      0  0]
    % R1 is an r*r upper triangular regular matrix
    % R2 is an r*(c-r) matrix where c is the number of standard parameters
    R1 = R(1:r,1:r);
    R2 = R(1:r,r+1:end);
    beta = R1\R2; % The zero rows of K correspond to independent columns of WP
    beta(abs(beta)<sqrt(eps)) = 0; % Eliminate numerical errors
    % W2 = W1*beta

    % Ensure the relationship holds
    W1 = W*E(:,1:r);
    W2 = W*E(:,r+1:end);
    assert(norm(W2 - W1*beta) < 1e-6,... 
            'Found relationship between W1 and W2 is not correct\n');
    
    % -----------------------------------------------------------------------
    % Minimal set of inertial parameters
    % -----------------------------------------------------------------------
    pi1 = E(:,1:r)'*rbt_df.link_params'; % Independent parameters
    pi2 = E(:,r+1:end)'*rbt_df.link_params'; % Dependent parameters
    
    base_params = pi1 + beta*pi2;



end