function base = DynamicBaseParamCalc(rbt_df, dyn)
    tic;
    disp("Calculating base parameter...");
    [r, P_X, P, Kd] = find_dyn_param_deps(rbt_df, dyn);
    disp('calc_base_param_: time ' + string(toc) + 'sec');
    
    n = length(rbt_df.bary_params);
    s = n-r;
    KG = [eye(r), Kd; zeros(s, r), eye(s)];
    eye_param_num = eye(n);
    inertia2beta = KG * eye_param_num(:, P)';

    prime_param = rbt_df.bary_params';
    beta_param = inertia2beta * prime_param;
    base_param = beta_param(1:r);


    P_b = P(1:r);
    H_b = dyn.H(:, P_b);

    base.base_num = r;
    base.prime_num = n;
    base.redundant_num = s;
    base.prime_param = prime_param;
    base.beta_param = beta_param;
    base.base_param = base_param(1:r);
    base.prime2beta = inertia2beta;
    base.H_b = H_b;
    base.H = dyn.H;
    
    
    
end


function [r, P_X, P, Kd] = find_dyn_param_deps(rbt_df, dyn)
    dof = rbt_df.dof;
    param_num = length(rbt_df.bary_params);
    sample_num = param_num*2;
    Z = zeros(dof * sample_num, param_num);
    
    funcHandle = matlabFunction(vpa(dyn.H),'Vars', {rbt_df.coordinates', rbt_df.d_coordinates', rbt_df.dd_coordinates'});
    
    for i = 1:sample_num
        q = rand([dof, 1])*2*pi - pi;
        dq = rand([dof, 1])*2*pi - pi;
        ddq = rand([dof, 1])*2*pi - pi;
    %     Z(((i-1)*dof+1):(i*dof), :) = subs(regressor, [rbt_df.coordinates, rbt_df.d_coordinates, rbt_df.dd_coordinates], [q,dq,ddq]);
        Z(((i-1)*dof+1):(i*dof), :) = funcHandle(q,dq,ddq);
    end
    %% For details, check the note N22-003 section 3
    r = rank(Z);
    [~, ~, P] = qr(Z, 'vector');
    
    [~, R] = qr(Z(:, P));
    R1 = R(1:r, 1:r);
    R2 = R(1:r, (r+1):end);
    eye_param_num = eye(param_num);
    
    P_X = [eye(r), R1\R2]*eye_param_num(:, P)'; %P_x mapping the parameter to the base set
    
    Kd = R1\R2;

end