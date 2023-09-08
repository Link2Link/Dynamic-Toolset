function geom = GeometryCalculation(rbt_df, static_model)
    if nargin > 0
        if nargin == 1
            static_model = false;
        end
        geom.static_model = static_model;
        
        tic;
        geom = calc_geom_(rbt_df, geom);
        geom = calc_functions_(rbt_df, geom);
        disp('GeometryCalculation: time ' + string(toc) + 'sec');
    end
end

function geom = calc_geom_(rbt_df, geom)
% calc forward knematic

    geom.T_0n = num2cell(zeros(1, rbt_df.frame_num));      % Tsb of each link
    geom.R = num2cell(zeros(1, rbt_df.frame_num));         % Tsb(1:3, 1:3)
    geom.p_n = num2cell(zeros(1, rbt_df.frame_num));       % Tsb(1:3, 4)
    geom.T_0nc = num2cell(zeros(1, rbt_df.frame_num));     % Ts_COM
    geom.p_c = num2cell(zeros(1, rbt_df.frame_num));       % COM positon
    
    geom.v_cw = num2cell(zeros(1, rbt_df.frame_num));      % speed of COM
    geom.w_b = num2cell(zeros(1, rbt_df.frame_num));
    
    Tsb = num2cell(zeros(1, rbt_df.frame_num));

    syms t real

    f = utils;
    if geom.static_model
        disp('static_model = true, Overlook the dynamic part');
    else
        disp('(default) static_model = false');
    end
    for num = 1:rbt_df.frame_num
        disp('Kinematic for Frame : ' + string(rbt_df.link_nums{num}));
        if num == 1
            geom.T_0n{num} = rbt_df.M{num};
            Tsb{num} = sym(eye(4));
            continue
        end

        Tsb{num} = Tsb{rbt_df.prev_link_num{num}+1} * f.matrix_exp_6(rbt_df.screw{num}, rbt_df.theta{num});
        geom.T_0n{num} = simplify(expand(Tsb{num}*rbt_df.M{num}));
        geom.R{num} = geom.T_0n{num}(1:3, 1:3);
        geom.p_n{num} = geom.T_0n{num}(1:3, 4); % joint position
        
        
        geom.T_0nc{num} = simplify(expand(geom.T_0n{num} * f.translation_transfmat(rbt_df.r_by_ml{num})));
        geom.p_c{num} = geom.T_0nc{num}(1:3, 4); % CoM position
        
        if ~geom.static_model   % if static_model, pass the differential kinematics
            v_cw_temp = diff(subs(geom.p_c{num}, rbt_df.subs_q2qt(1,:), rbt_df.subs_q2qt(2,:)), t);
            v_cw_temp = subs(v_cw_temp, [rbt_df.subs_dqt2dq(1,:), rbt_df.subs_qt2q(1,:)], [rbt_df.subs_dqt2dq(2,:), rbt_df.subs_qt2q(2,:)]);
            geom.v_cw{num} = simplify(expand(v_cw_temp));

            R_t = subs(geom.R{num}, rbt_df.subs_q2qt(1,:), rbt_df.subs_q2qt(2,:));
            dR_t = diff(R_t, t);
            dR = subs(dR_t, [rbt_df.subs_dqt2dq(1,:), rbt_df.subs_qt2q(1,:)], [rbt_df.subs_dqt2dq(2,:), rbt_df.subs_qt2q(2,:)]);

            geom.w_b{num} = simplify(expand(f.so32vec((geom.R{num})'*dR)));
        end

    end
end



function geom = calc_functions_(rbt_df, geom)
    geom.p_n_func = num2cell(zeros(1, rbt_df.frame_num));
    for num = 1:rbt_df.frame_num
        geom.p_n_func{num} = matlabFunction(geom.p_n{num}, 'Vars', rbt_df.coordinates);
    end
end