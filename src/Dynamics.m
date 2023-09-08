function dyn = Dynamics(rbt_df, geom, gravity_vec)
    dyn = struct();
    gravity_vec = gravity_vec(:)';
    dyn.gravity_vec    = gravity_vec;
    tic;
    if geom.static_model
        disp('static_model = true, overlook the dynamic effect.');
    else
        disp('(default), static_model = false');
    end
    disp("Calculating Lagrangian...")

    p_e = 0;
    k_e = 0;

    for num = cell2mat(rbt_df.link_nums(2:end)')
        k_e_n = 0;
        if rbt_df.use_inertia{num+1}
            disp("Calculating the link kinetic energy of " + num + "/" + rbt_df.link_nums{end});
            p_e = p_e - rbt_df.m{num+1} * gravity_vec * geom.p_c{num+1};


            if ~geom.static_model
                k_e_n = rbt_df.m{num+1} * geom.v_cw{num+1}' * (geom.v_cw{num+1}) / 2 + ...
                geom.w_b{num+1}' * rbt_df.I_by_Llm{num+1} * geom.w_b{num+1} / 2;

                k_e_n = prod(factor(expand(k_e_n) - subs(expand(k_e_n * rbt_df.m{num+1}), rbt_df.m{num+1}, 0)/rbt_df.m{num+1}), "All");
            end

        end
        k_e = k_e+ k_e_n; 
    end
    dyn.Lagrange.K = k_e;
    dyn.Lagrange.P = p_e;
    % Lagrangian
    L = k_e - p_e;
    dyn.Lagrange.L = L;
    dyn.tau = [];

    disp("Calculating joint torques...");
    syms t real

    for content = [rbt_df.coordinates; rbt_df.d_coordinates]
        q = content(1);
        dq = content(2);
        disp("tau of " + string(q));
        dk_ddq = diff(k_e, dq);
        dk_ddq_t = subs(dk_ddq, [rbt_df.subs_q2qt(1,:), rbt_df.subs_dq2dqt(1,:)], [rbt_df.subs_q2qt(2,:), rbt_df.subs_dq2dqt(2,:)]);
        dk_ddq_dt_t = diff(dk_ddq_t, t);

        dk_ddq_dt = subs(dk_ddq_dt_t, [rbt_df.subs_ddqt2ddq(1,:), rbt_df.subs_dqt2dq(1,:), rbt_df.subs_qt2q(1,:)], [rbt_df.subs_ddqt2ddq(2,:), rbt_df.subs_dqt2dq(2,:), rbt_df.subs_qt2q(2,:)]);

        dL_dq = diff(L, q);

        dyn.tau = [dyn.tau, expand(dk_ddq_dt - dL_dq)];
    end
    dyn.tau = dyn.tau';
    dyn.tau_MCG = dyn.tau;
    
    disp("Adding springs...");
    dyn.tau_spring = sym(zeros(length(rbt_df.d_coordinates), 1));
    for i = 1:rbt_df.frame_num
        q = rbt_df.q_for_frame(i);
        dq = rbt_df.dq_for_frame(i);

%         if rbt_df.use_friction{i}          % not include friction
%             tau_f = sign(dq) * rbt_df.Fc{i} + dq * rbt_df.Fv{i};
%             for a = 1:length(rbt_df.d_coordinates)
%                 dq_da = diff(dq, rbt_df.d_coordinates(a));
%                 obj.tau(a) = obj.tau(a) + dq_da * tau_f;
%             end
%         end
        
        if ~isnan(rbt_df.spring_dl{i})
            tau_s = sum(rbt_df.spring_dl{i} .* rbt_df.spring_param{i});
            for a = 1:length(rbt_df.d_coordinates)
                dq_da = diff(dq, rbt_df.d_coordinates(a));
                dyn.tau_spring(a) = dyn.tau_spring(a) + dq_da * tau_s;
            end
        end
    end
    dyn.tau = dyn.tau_MCG + dyn.tau_spring;
      
%     disp('final simplifying');        % consume too many time, skip it
%     dyn.tau = simplify(dyn.tau);
    
    disp('calc_dyn_: time ' + string(toc) + 'sec');
    dyn.H = Dynamics_calc_regressor(rbt_df, dyn);
    dyn.inertial_param = rbt_df.bary_params';
end


function H = Dynamics_calc_regressor(rbt_df, dyn)
    tic;
    disp("Calculating regressor...");
    [H, ~] = equationsToMatrix(dyn.tau, rbt_df.bary_params);
    disp('calc_regressor_: time ' + string(toc) + 'sec');
end