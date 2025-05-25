function geom = GeometryCalculation(rbt_df, static_model)
%GeometryCalculation 计算运动学
%   rbt_df : DefineRobot的结果
%   static_model : 是否只考虑重力，不考虑M和C
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

    geom.T_0n = num2cell(zeros(1, rbt_df.frame_num));      % 末端位姿
    geom.R = num2cell(zeros(1, rbt_df.frame_num));         % 末端姿态
    geom.p_n = num2cell(zeros(1, rbt_df.frame_num));       % 末端位置
    geom.T_0nc = num2cell(zeros(1, rbt_df.frame_num));     % COM系位姿矩阵
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

        % 指数映射的乘积
        Tsb{num} = Tsb{rbt_df.prev_link_num{num}+1} * f.matrix_exp_6(rbt_df.screw{num}, rbt_df.theta{num});
        % POE公式
        geom.T_0n{num} = simplify(expand(Tsb{num}*rbt_df.M{num}));
        % 各个坐标系的姿态
        geom.R{num} = geom.T_0n{num}(1:3, 1:3);
        % 各个坐标系的位置
        geom.p_n{num} = geom.T_0n{num}(1:3, 4); 
        % COM系的位姿
        geom.T_0nc{num} = simplify(expand(geom.T_0n{num} * f.translation_transfmat(rbt_df.r_by_ml{num})));
        % COM点的位置
        geom.p_c{num} = geom.T_0nc{num}(1:3, 4); % CoM position
        
        if ~geom.static_model   % if static_model, pass the differential kinematics
            % 质心速度（表达在基坐标系）
            v_cw_temp = diff(subs(geom.p_c{num}, rbt_df.subs_q2qt(1,:), rbt_df.subs_q2qt(2,:)), t);
            v_cw_temp = subs(v_cw_temp, [rbt_df.subs_dqt2dq(1,:), rbt_df.subs_qt2q(1,:)], [rbt_df.subs_dqt2dq(2,:), rbt_df.subs_qt2q(2,:)]);
            geom.v_cw{num} = simplify(expand(v_cw_temp));   

            % 姿态矩阵的微分
            R_t = subs(geom.R{num}, rbt_df.subs_q2qt(1,:), rbt_df.subs_q2qt(2,:));
            dR_t = diff(R_t, t);
            dR = subs(dR_t, [rbt_df.subs_dqt2dq(1,:), rbt_df.subs_qt2q(1,:)], [rbt_df.subs_dqt2dq(2,:), rbt_df.subs_qt2q(2,:)]);
            
            % 体坐标旋量wb so3
            geom.w_b{num} = simplify(expand(f.so32vec((geom.R{num})'*dR)));
        end

    end
end



function geom = calc_functions_(rbt_df, geom)
    geom.p_n_func = num2cell(zeros(1, rbt_df. frame_num));
    for num = 1:rbt_df.frame_num
        geom.p_n_func{num} = matlabFunction(geom.p_n{num}, 'Vars', rbt_df.coordinates);
    end
end