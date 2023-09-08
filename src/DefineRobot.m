function rbt = DefineRobot(model_name,params)
%DEFINEROBOT Define the robot from model file
%   model_name : name of the robot 
%   params : config data
rbt = struct();
rbt.name = model_name;
rbt.frame_num = height(params);
rbt.link_nums = params(:, 1);           % link number
rbt.prev_link_num = params(:, 2);       % previous link
rbt.succ_link_num = params(:, 3);       % successor link
rbt.shat = params(:, 4);                % screw axis unit direction
rbt.rs = params(:, 5);                  % rs
rbt.M_R = params(:, 6);                 % rotation matrix of M
rbt.theta = params(:, 7);               % screw angle
rbt.h = params(:, 8);                   % screw pitch
rbt.use_inertia = params(:, 9);         % link inertia
rbt.spring_dl = params(:, 10);          % spring

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
        screw{num} = [ws{num}; vs{num}];
        M = [M, {[rbt.M_R{num}, rbt.rs{num}; [0,0,0,1]]}];

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
    rbt.m = num2cell(zeros(1, rbt.frame_num));              % mass            
    rbt.l = num2cell(zeros(1, rbt.frame_num));              % COM in link frame
    rbt.r = num2cell(zeros(1, rbt.frame_num));              % link origin in COM frame
    rbt.r_by_ml = num2cell(zeros(1, rbt.frame_num));        % r represent by L/m
    rbt.L_vec = num2cell(zeros(1, rbt.frame_num));          % inertial tenser vecotr in link frame
    rbt.I_vec = num2cell(zeros(1, rbt.frame_num));          % inertial tenser vector in COM frame
    rbt.L_mat = num2cell(zeros(1, rbt.frame_num));          
    rbt.I_mat = num2cell(zeros(1, rbt.frame_num));
    rbt.I_by_Llm = num2cell(zeros(1, rbt.frame_num));
    rbt.spring_formula = num2cell(zeros(1, rbt.frame_num));
    rbt.spring_param = num2cell(zeros(1, rbt.frame_num));
    rbt.std_params = [];                                    % inertial parameter in link frame
    rbt.bary_params = [];                                   % inertial parameter in COM frame
    
    spring_num = 0;
    for num = 1:(length(rbt.link_nums)-1)
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
            % inertial parameter in link frame
            rbt.bary_params = [rbt.bary_params, rbt.L_vec{num+1}];
            rbt.bary_params = [rbt.bary_params, rbt.l{num+1}];
            rbt.bary_params = [rbt.bary_params, [rbt.m{num+1}]];
            
            % inertial parameter in mass center frame
            rbt.std_params = [rbt.std_params, rbt.I_vec{num+1}];
            rbt.std_params = [rbt.std_params, rbt.r{num+1}];
            rbt.std_params = [rbt.std_params, [rbt.m{num+1}]];
        end
    end

    rbt.inertial_param_number = length(rbt.bary_params);
    for num = 1:(length(rbt.link_nums)-1)
        if ~isnan(rbt.spring_dl{num+1})
            rbt.bary_params = [rbt.bary_params, [rbt.spring_param{num+1}]];
            rbt.std_params = [rbt.std_params, [rbt.spring_param{num+1}]];
        end
    end
    rbt.spring_param_number = spring_num;

end


function rbt = gen_coordinates_(rbt)
rbt.coordinates = [];
rbt.coordinates_joint_type = [];

% number of frames
for num = 1:length(rbt.link_nums)
    for s = symvar(rbt.theta{num})
        if ~ismember(s, rbt.coordinates)
            rbt.coordinates = [rbt.coordinates, s];
            rbt.coordinates_joint_type = [rbt.coordinates_joint_type, rbt.joint_type(num)];
        end
    end
end
rbt.dof = length(rbt.coordinates);

rbt.d_coordinates = [];
rbt.dd_coordinates = [];
rbt.coordinates_t = [];
syms t real
for co = rbt.coordinates
    rbt.d_coordinates = [rbt.d_coordinates, sym("d"+string(co), 'real')];
    rbt.dd_coordinates = [rbt.dd_coordinates, sym("dd"+string(co), 'real')];

    syms(string(co)+"t(t)");
    rbt.coordinates_t = [rbt.coordinates_t, eval(string(co)+"t(t)")];
end

rbt.d_coordinates_t = [];
rbt.dd_coordinates_t = [];
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
    %             q = NaN;
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

