function rbt = DefineRobot(model_name,params)
%DEFINEROBOT 根据给定配置生成所需基础符号
%   model_name : 模型名
%   params : 模型参数表
rbt = struct();                         % 基础结构体，后续所有中间变量和最终结果都保存在此结构体中
rbt.name = model_name;                  % 模型名
rbt.frame_num = height(params);         % 模型中坐标系数量
rbt.link_nums = params(:, 1);           % 连杆编号
rbt.prev_link_num = params(:, 2);       % 前一级连杆编号
rbt.succ_link_num = params(:, 3);       % 后一级连杆编号
rbt.shat = params(:, 4);                % 旋量轴朝向(单位向量)
rbt.rs = params(:, 5);                  % 旋量轴上一点(基坐标系下表示)
rbt.M_R = params(:, 6);                 % 初始位姿下坐标系的姿态矩阵(基坐标系下表示)
rbt.theta = params(:, 7);               % 关节自变量
rbt.h = params(:, 8);                   % 节距(纯旋转h=0, 纯平移h=inf, 螺旋运动h为有限值)
rbt.use_inertia = params(:, 9);         % 是否考虑本连杆的动力学效应
rbt.spring_dl = params(:, 10);          % 关节弹簧模型

rbt = DefineRobot_gentransfm(rbt);
rbt = DefineRobot_genparams(rbt);
rbt = gen_coordinates_(rbt);
end

function rbt = DefineRobot_gentransfm(rbt)
    M = {};
    ws = num2cell(zeros(3, rbt.frame_num),1);
    vs = num2cell(zeros(3, rbt.frame_num),1);
    screw = num2cell(zeros(6, rbt.frame_num),1);
    rbt.joint_type = [];
    
    for num = 1:length(rbt.link_nums)
        if isinf(rbt.h{num})
            ws{num} = [0;0;0];
            h_temp = 1;
        else
            ws{num} = rbt.shat{num};
            h_temp = rbt.h{num};
        end
        vs{num} = -cross(ws{num}, rbt.rs{num}) + h_temp*rbt.shat{num};
        screw{num} = [ws{num}; vs{num}];                                    % 计算旋量
        M = [M, {[rbt.M_R{num}, rbt.rs{num}; [0,0,0,1]]}];

        % 根据旋量判断关节类型
        % TODO: 此方法仅适用于1自由度低幅，下一大版本考虑改掉这里以支持2自由度和自由度低副
        if norm(screw{num})<eps
            rbt.joint_type = [rbt.joint_type, "F"];  % Fixed
        elseif isinf(rbt.h{num})
            rbt.joint_type = [rbt.joint_type, "P"];  % Prismatic
        elseif rbt.h{num}==0
           rbt.joint_type = [rbt.joint_type, "R"];  % Revolute
        else
           rbt.joint_type = [rbt.joint_type, "A"];  % Assitive
        end
    
    end
    rbt.M = M;
    rbt.ws = ws;
    rbt.vs = vs;
    rbt.screw = screw;
end



function rbt = DefineRobot_genparams(rbt)

    % inertial parameter
    rbt.m = num2cell(zeros(1, rbt.frame_num));              % 连杆质量            
    rbt.l = num2cell(zeros(1, rbt.frame_num));              % 质量矩
    rbt.r = num2cell(zeros(1, rbt.frame_num));              % 关节坐标系下的COM位置
    rbt.r_by_ml = num2cell(zeros(1, rbt.frame_num));        % r 表示为 l/m
    rbt.L_vec = num2cell(zeros(1, rbt.frame_num));          % 关节系下的惯性张量(向量)
    rbt.I_vec = num2cell(zeros(1, rbt.frame_num));          % COM系下的惯性张量(向量)
    rbt.L_mat = num2cell(zeros(1, rbt.frame_num));          % 关节系下的惯性张量(矩阵)
    rbt.I_mat = num2cell(zeros(1, rbt.frame_num));          % COM系下的惯性张量(矩阵)
    rbt.I_by_Llm = num2cell(zeros(1, rbt.frame_num));       % 使用L、l、m表示I
    rbt.spring_formula = num2cell(zeros(1, rbt.frame_num)); % 弹簧模型
    rbt.spring_param = num2cell(zeros(1, rbt.frame_num));   % 弹簧参数
    rbt.std_params = [];                                    % COM系下的惯性参数
    rbt.bary_params = [];                                   % 关节系下的惯性参数
    
    spring_num = 0;
    for num = 1:(length(rbt.link_nums)-1)
        % 使用num+1是跳过第一个坐标系，第一个坐标系默认为基座坐标
        % TODO:这里的基础假设是固定基，下一版本考虑支持浮动基
        rbt.m{num+1} = sym('m'+string(num), 'real');
        rbt.l{num+1} = [sym('l'+string(num)+'x', 'real'), sym('l'+string(num)+'y', 'real'), sym('l'+string(num)+'z', 'real')];
        rbt.r{num+1} = [sym('r'+string(num)+'x', 'real'), sym('r'+string(num)+'y', 'real'), sym('r'+string(num)+'z', 'real')];
        rbt.I_vec{num+1} = [sym('I'+string(num)+'xx', 'real'), sym('I'+string(num)+'xy', 'real'), sym('I'+string(num)+'xz', 'real'), sym('I'+string(num)+'yy', 'real'), sym('I'+string(num)+'yz', 'real'), sym('I'+string(num)+'zz', 'real')];
        rbt.L_vec{num+1} = [sym('L'+string(num)+'xx', 'real'), sym('L'+string(num)+'xy', 'real'), sym('L'+string(num)+'xz', 'real'), sym('L'+string(num)+'yy', 'real'), sym('L'+string(num)+'yz', 'real'), sym('L'+string(num)+'zz', 'real')];

        f = utils;
        rbt.I_mat{num+1} = f.inertia_vec2tensor(rbt.I_vec{num+1});
        rbt.L_mat{num+1} = f.inertia_vec2tensor(rbt.L_vec{num+1});

        rbt.r_by_ml{num+1} = f.ml2r(rbt.m{num+1}, rbt.l{num+1});
        rbt.I_by_Llm{num+1} = f.Lmr2I(rbt.L_mat{num+1}, rbt.m{num+1}, rbt.r_by_ml{num+1});

        if ~isnan(rbt.spring_dl{num+1})
           spring_param_number = width(rbt.spring_dl{num+1});
           rbt.spring_param{num+1} = sym('K' + string(num), [1,spring_param_number], 'real');
           rbt.spring_formula{num+1} = sum(rbt.spring_dl{num+1} .* rbt.spring_param{num+1});
           spring_num = spring_num + spring_param_number;
        end
    end
    

    
    for num = 1:(length(rbt.link_nums)-1)
        if rbt.use_inertia{num+1}
            % 关节坐标系惯性参数
            rbt.bary_params = [rbt.bary_params, rbt.L_vec{num+1}];
            rbt.bary_params = [rbt.bary_params, rbt.l{num+1}];
            rbt.bary_params = [rbt.bary_params, [rbt.m{num+1}]];
            
            % COM系惯性参数
            rbt.std_params = [rbt.std_params, rbt.I_vec{num+1}];
            rbt.std_params = [rbt.std_params, rbt.r{num+1}];
            rbt.std_params = [rbt.std_params, [rbt.m{num+1}]];
        end
    end

    % 惯性参数数量等于10x连杆数量
    rbt.inertial_param_number = length(rbt.bary_params);

    % 在惯性参数后面添加额外的弹簧补偿参数
    for num = 1:(length(rbt.link_nums)-1)
        if ~isnan(rbt.spring_dl{num+1})
            rbt.bary_params = [rbt.bary_params, [rbt.spring_param{num+1}]];
            rbt.std_params = [rbt.std_params, [rbt.spring_param{num+1}]];
        end
    end
    % 弹簧参数数量根据config中模型给定
    rbt.spring_param_number = spring_num;

end


function rbt = gen_coordinates_(rbt)
rbt.coordinates = [];
rbt.coordinates_joint_type = [];

% number of frames
for num = 1:length(rbt.link_nums)
    % 考虑一个关节有多个自变量的情况
    for s = symvar(rbt.theta{num})
        if ~ismember(s, rbt.coordinates)
            rbt.coordinates = [rbt.coordinates, s];
            rbt.coordinates_joint_type = [rbt.coordinates_joint_type, rbt.joint_type(num)];
        end
    end
end
rbt.dof = length(rbt.coordinates);

rbt.d_coordinates = [];     %存放速度符号dq
rbt.dd_coordinates = [];    %存放加速度符号ddq
rbt.coordinates_t = [];     %存放带时间符号q(t)
syms t real
for co = rbt.coordinates
    rbt.d_coordinates = [rbt.d_coordinates, sym("d"+string(co), 'real')];
    rbt.dd_coordinates = [rbt.dd_coordinates, sym("dd"+string(co), 'real')];

    syms(string(co)+"t(t)");
    rbt.coordinates_t = [rbt.coordinates_t, eval(string(co)+"t(t)")];
end

rbt.d_coordinates_t = [];   %存放q(t)对时间的微分
rbt.dd_coordinates_t = [];  %存放q(t)对时间的二次微分
for co_t = rbt.coordinates_t
    rbt.d_coordinates_t = [rbt.d_coordinates_t, diff(co_t, t)];
    rbt.dd_coordinates_t = [rbt.dd_coordinates_t, diff(rbt.d_coordinates_t(end), t)];
end

rbt.subs_q2qt = [];
rbt.subs_dq2dqt = [];
rbt.subs_ddq2ddqt = [];
rbt.subs_qt2q = [];
rbt.subs_dqt2dq = [];
rbt.subs_ddqt2ddq = [];


% 构造q、dq、ddq和q(t)、dq(t)、ddq(t)之间的互换对应
rbt.subs_q2qt = [rbt.coordinates; rbt.coordinates_t];
rbt.subs_dq2dqt = [rbt.d_coordinates; rbt.d_coordinates_t];
rbt.subs_ddq2ddqt = [rbt.dd_coordinates; rbt.dd_coordinates_t];
rbt.subs_qt2q = [rbt.coordinates_t; rbt.coordinates];
rbt.subs_dqt2dq = [rbt.d_coordinates_t; rbt.d_coordinates];
rbt.subs_ddqt2ddq = [rbt.dd_coordinates_t; rbt.dd_coordinates];

rbt.q_for_frame = sym(zeros(1, rbt.frame_num));
rbt.dq_for_frame = sym(zeros(1, rbt.frame_num));
rbt.ddq_for_frame = sym(zeros(1, rbt.frame_num));

for i = 1:rbt.frame_num
    % q = NaN;
    if rbt.joint_type(i) == "P" || rbt.joint_type(i) == "R" || rbt.joint_type(i) == "A"
        q = rbt.theta{i};
    else
        continue
    end
    rbt.q_for_frame(i) = q;
    qt = subs(q, rbt.subs_q2qt(1,:), rbt.subs_q2qt(2,:)); % q2qt
    dqt = diff(qt, t);
    dq = subs(dqt, rbt.subs_dqt2dq(1,:), rbt.subs_dqt2dq(2,:)); % dqt2dq

    rbt.dq_for_frame(i) = dq;

    ddqt = diff(dqt, t);
    ddq = subs(ddqt, rbt.subs_ddqt2ddq(1,:), rbt.subs_ddqt2ddq(2,:)); % ddqt2ddq

    rbt.ddq_for_frame(i) = ddq;
end

end

