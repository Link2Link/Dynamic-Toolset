function base = DynamicBaseParamCalc(rbt_df, dyn)
    tic;
    disp("Calculating base parameter...");
    [r, P_X, P, Kd] = find_dyn_param_deps(rbt_df, dyn);
    disp('calc_base_param_: time ' + string(toc) + 'sec');
    
    n = length(rbt_df.bary_params);
    s = n-r;
    KG = [eye(r), Kd; zeros(s, r), eye(s)];
    eye_param_num = eye(n);
    inertia2beta = KG * eye_param_num(:, P)'; % 原始惯性参数到beta集合的映射  详细原来查看N22-003 section 5

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
    base.H_b = collect(expand(simplify(H_b)));
    base.H = dyn.H;
    
    
    
end


function [r, P_X, P, Kd] = find_dyn_param_deps(rbt_df, dyn)
    dof = rbt_df.dof;
    param_num = length(rbt_df.bary_params);
    % 采样数量等于参数数量的两倍，确保能够拿到极大线性无关组
    sample_num = param_num*2;
    Z = zeros(dof * sample_num, param_num);
    H = vpa(dyn.H);
    if isa(dyn.gravity_vec, 'sym')
        H = subs(H, dyn.gravity_vec, rand(1,3)); % 将重力向量换为数值再进行惯性参数分解
    end
    funcHandle = matlabFunction(H,'Vars', {rbt_df.coordinates', rbt_df.d_coordinates', rbt_df.dd_coordinates'});
    % 固定随机数种子，确保每次生成结果相同
    rng(123);   
    for i = 1:sample_num
        q = rand([dof, 1])*2*pi - pi;
        dq = rand([dof, 1])*2*pi - pi;
        ddq = rand([dof, 1])*2*pi - pi;
    %     Z(((i-1)*dof+1):(i*dof), :) = subs(regressor, [rbt_df.coordinates, rbt_df.d_coordinates, rbt_df.dd_coordinates], [q,dq,ddq]);
        Z(((i-1)*dof+1):(i*dof), :) = funcHandle(q,dq,ddq);
    end
    %% 使用QR分解计算最小惯性参数集, 详细原理查看 N22-003 section 5
    r = rank(Z);
    [~, ~, P] = qr(Z, 'vector');
    
    [~, R] = qr(Z(:, P));
    R1 = R(1:r, 1:r);
    R2 = R(1:r, (r+1):end);
    eye_param_num = eye(param_num);
    
    P_X = [eye(r), R1\R2]*eye_param_num(:, P)'; %P_x mapping the parameter to the base set
    Kd = R1\R2;
    
    % 微小项强制给0
    P_X(abs(P_X) < 1e-10) = 0;
    Kd(abs(Kd) < 1e-10) = 0;

end