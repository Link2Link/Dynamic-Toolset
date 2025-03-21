function rbt_df = DefineRobot_MDH(model_name,params)
%DEFINEROBOT Define the robot from model file
%   model_name : name of the robot 
%   mode : mode of Kinematics
%   params : config data
    tic;
    rbt_df = struct();
    rbt_df.name = model_name;
    rbt_df.frame_num = height(params);
    rbt_df.dof = rbt_df.frame_num-1;
    rbt_df.link_nums = params(:, 1);           % link number
    rbt_df.prev_link_num = params(:, 2);       % previous link
    rbt_df.succ_link_num = params(:, 3);       % successor link
    rbt_df.a = params(:, 4);                   % a of mdh
    rbt_df.alpha = params(:, 5);               % alpha of mdh
    rbt_df.d = params(:, 6);                   % d of mdh
    rbt_df.theta = params(:, 7);               % theta of mdh
    rbt_df.offset = params(:, 8);              % theta offset
    rbt_df.joint_type = params(:, 9);          % Joint Type
    rbt_df.use_inertia = params(:, 10);        % link inertia
    
    
    rbt_df = DefineRobot_genParams_MDH(rbt_df);

end

function rbt_df = DefineRobot_genParams_MDH(rbt_df)

    % inertial parameter
    rbt_df.m = num2cell(zeros(1, rbt_df.dof));              % mass
    rbt_df.r = num2cell(zeros(1, rbt_df.dof));              % COM in link frame
    rbt_df.Ic_vac = num2cell(zeros(1, rbt_df.dof));         % inertial tenser vecotr in COM frame
    rbt_df.l = num2cell(zeros(1, rbt_df.dof));              % first moment of inertial
    rbt_df.L_vec = num2cell(zeros(1, rbt_df.dof));          % inertial tenser vecotr in link frame
    rbt_df.L_mat = num2cell(zeros(1, rbt_df.dof));          % inertial tenser matrix in link frame
    rbt_df.link_params = [];                             % inertial parameter in link frame
    rbt_df.std_params = [];                              % standard inertial parameters
    rbt_df.coordinates = sym(zeros(1, rbt_df.dof));               
    rbt_df.d_coordinates = sym(zeros(1, rbt_df.dof));
    rbt_df.dd_coordinates = sym(zeros(1, rbt_df.dof));
    
    for num = 1:rbt_df.dof
        rbt_df.m{num} = sym('m'+string(num), 'real');
        rbt_df.r{num} = [sym('r'+string(num)+'x', 'real'), sym('r'+string(num)+'y', 'real'), sym('r'+string(num)+'z', 'real')];
        rbt_df.Ic_vec{num} = [sym('Ic'+string(num)+'xx', 'real'), sym('Ic'+string(num)+'xy', 'real'), sym('Ic'+string(num)+'xz', 'real'), sym('Ic'+string(num)+'yy', 'real'), sym('Ic'+string(num)+'yz', 'real'), sym('Ic'+string(num)+'zz', 'real')];
        rbt_df.l{num} = [sym('l'+string(num)+'x', 'real'), sym('l'+string(num)+'y', 'real'), sym('l'+string(num)+'z', 'real')];
        rbt_df.L_vec{num} = [sym('L'+string(num)+'xx', 'real'), sym('L'+string(num)+'xy', 'real'), sym('L'+string(num)+'xz', 'real'), sym('L'+string(num)+'yy', 'real'), sym('L'+string(num)+'yz', 'real'), sym('L'+string(num)+'zz', 'real')];

        f = utils;
        % the first moment of inertial in urdf is Ic , no need minus
        rbt_df.Ic_mat{num} = f.inertia_vec2tensor(rbt_df.Ic_vec{num}); % inertial vector to matrix
        rbt_df.L_mat{num} = f.inertia_vec2tensor(rbt_df.L_vec{num});   % inertial vector to matrix 
        rbt_df.link_params = [rbt_df.link_params, rbt_df.L_vec{num}, rbt_df.l{num}, rbt_df.m{num}]; % [Lxx Lxy Lxz Lyy Lyz Lzz lx ly lz m]
        rbt_df.std_params = [rbt_df.std_params, rbt_df.Ic_vec{num}, rbt_df.r{num}, rbt_df.m{num}]; % [Icxx Icxy Icxz Icyy Icyz Iczz rx ry rz m]
    end

    rbt_df.inertial_param_number = length(rbt_df.std_params);
    
    for i = 1:rbt_df.dof
        if rbt_df.joint_type{i+1} ~= 3
            rbt_df.coordinates(i) = rbt_df.theta{i+1};    
            rbt_df.d_coordinates(i) = sym("d"+string(rbt_df.theta{i+1}), 'real');
            rbt_df.dd_coordinates(i) = sym("dd"+string(rbt_df.theta{i+1}), 'real');
        end
        
    end
    
    disp('link_params:'+string(size(rbt_df.link_params,2)) + ' ,DefineRobot: time ' + string(toc) + 'sec');
end
