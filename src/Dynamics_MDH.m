function dyn = Dynamics_MDH(rbt_df, geom, gravity_vec)
    tic;
    format long;
    Rn_n1 = geom.Rn_n1;
    Pn_n1 = geom.Pn_n1;
    dq = rbt_df.d_coordinates;
    ddq = rbt_df.dd_coordinates;
    w = cell(1,rbt_df.frame_num); % Angular velocity
    dw = cell(1,rbt_df.frame_num); % Angular acceleration
    dv = cell(1,rbt_df.frame_num); % Linear acceleration
    % Centroid inertia modeling
    dv_c_ir = cell(1,rbt_df.frame_num); % Centroid linear acceleration
    Fc_ir = cell(1,rbt_df.frame_num); % Centroid force
    Nc_ir = cell(1,rbt_df.frame_num); % Centroid Torque
    f_ir = cell(1,rbt_df.frame_num); % Motor output torque
    tau_ir = cell(1,rbt_df.frame_num); % Motor output force
    motor_t_ir = sym(zeros(rbt_df.dof,1)); % Final motor force/torque
    
    % Joint inertia modeling
    Fc_l = cell(1,rbt_df.frame_num); % Centroid Force
    f_l = cell(1,rbt_df.frame_num); % Motor output torque
    tau_l = cell(1,rbt_df.frame_num); % Motor output force
    motor_t_l = sym(zeros(rbt_df.dof,1)); % Final motor force/torque
    
    % Initialization
    w{1} = zeros(3,1);
    dw{1} = zeros(3,1);
    dv{1} = gravity_vec;  % Consider gravity
    Fc_ir{1} = zeros(3,1);
    Nc_ir{1} = zeros(3,1);

    Fc_l{1} = zeros(3,1);
    
    f_ir{rbt_df.frame_num} = zeros(3,1); 
    tau_ir{rbt_df.frame_num} = zeros(3,1);
    Pn_n1{rbt_df.frame_num} = zeros(3,1);
    Rn_n1{rbt_df.frame_num} = eye(3);

    f_l{rbt_df.frame_num} = zeros(3,1); 
    tau_l{rbt_df.frame_num} = zeros(3,1);
    %% Newton-Euler dynamics recursion -> Centroid inertia
    disp('COM inertial Newton-Euler')
    Zaxis = [0;0;1];% Rotary shaft
    % Forward pass
    for i=1:rbt_df.frame_num-1
        joint_type = rbt_df.joint_type{i+1};
        
        % Velocity/acceleration recursion
        if joint_type == 0  % Revolute joint
            w{i+1} = Rn_n1{i}' * w{i} + Zaxis*dq(i);
            dw{i+1} = Rn_n1{i}'*dw{i} + cross(Rn_n1{i}'*w{i}, Zaxis*dq(i)) + Zaxis*ddq(i);           
        else  % Prismatic joint
            w{i+1} = Rn_n1{i}' * w{i};
            dw{i+1} = Rn_n1{i}' * dw{i};
        end
        dv{i+1} = Rn_n1{i}'*(dv{i} + cross(dw{i}, Pn_n1{i}) + cross(w{i},cross(w{i},Pn_n1{i})));
        % Centroid acceleration
        dv_c_ir{i+1} = dv{i+1} + cross(dw{i+1}, rbt_df.r{i}') + cross(w{i+1}, cross(w{i+1}, rbt_df.r{i}'));
        % Inertial force
        Fc_ir{i+1} = rbt_df.m{i} * dv_c_ir{i+1};
        Nc_ir{i+1} = rbt_df.Ic_mat{i} * dw{i+1} + cross(w{i+1}, rbt_df.Ic_mat{i} * w{i+1});

        disp('Forward Pass:'+ string(i))
    end


    % Backward pass
    for i = rbt_df.frame_num:-1:2
        joint_type = rbt_df.joint_type{i}; 

        % Force/torque balance
        f_ir{i-1} = Rn_n1{i} * f_ir{i} + Fc_ir{i};
        tau_ir{i-1} = Rn_n1{i} * tau_ir{i} + cross(rbt_df.r{i-1}', Fc_ir{i}) + Nc_ir{i} + cross(Pn_n1{i},Rn_n1{i} * f_ir{i});

        % Extract joint torque
        if joint_type == 0
            motor_t_ir(i-1) = tau_ir{i-1}(3);  % Revolute joint takes z-axis torque
        else
            motor_t_ir(i-1) = f_ir{i-1}(3);    % Prismatic joint takes z-axis force
        end
        disp('Backward Pass:'+ string(i-1))
    end

    %% Newton-Euler dynamics recursion -> Joint inertia
    disp('Link inertial Newton-Euler')
    Zaxis = [0;0;1];% Rotary shaft
    % Forward pass
    for i=1:rbt_df.frame_num-1
        joint_type = rbt_df.joint_type{i+1};
        
        % Velocity/acceleration recursion
        if joint_type == 0  % Revolute joint
            w{i+1} = Rn_n1{i}' * w{i} + Zaxis*dq(i);
            dw{i+1} = Rn_n1{i}'*dw{i} + cross(Rn_n1{i}'*w{i}, Zaxis*dq(i)) + Zaxis*ddq(i);           
        else  % Prismatic joint
            w{i+1} = Rn_n1{i}' * w{i};
            dw{i+1} = Rn_n1{i}' * dw{i};
        end
        dv{i+1} = Rn_n1{i}'*(dv{i} + cross(dw{i}, Pn_n1{i}) + cross(w{i},cross(w{i},Pn_n1{i})));
        % Inertial force
        Fc_l{i+1} = rbt_df.m{i} * dv{i+1} + cross(dw{i+1},rbt_df.l{i}') + cross(w{i+1}, cross(w{i+1}, rbt_df.l{i}'));

        disp('Forward Pass:'+ string(i))
    end


    % Backward pass
    for i = rbt_df.frame_num:-1:2
        joint_type = rbt_df.joint_type{i}; 

        % Force/torque balance
        f_l{i-1} = Rn_n1{i} * f_l{i} + Fc_l{i};
        tau_l{i-1} = rbt_df.L_mat{i-1}*dw{i} + cross(w{i}, rbt_df.L_mat{i-1}*w{i}) + cross(rbt_df.l{i-1}',dv{i}) + Rn_n1{i} * tau_l{i} + cross(Pn_n1{i},Rn_n1{i} * f_l{i});

        % Extract joint torque
        if joint_type == 0
            motor_t_l(i-1) = tau_l{i-1}(3);  % Revolute joint takes z-axis torque
        else
            motor_t_l(i-1) = f_l{i-1}(3);    % Prismatic joint takes z-axis force
        end
        disp('Backward Pass:'+ string(i-1))
    end

    dyn.w = w; % Angular velocity
    dyn.dw = dw; % Angular acceleration
    dyn.dv = dv; % Linear acceleration
    dyn.dv_c_ir = dv_c_ir; % Linear acceleration
    dyn.Fc_ir = Fc_ir; % Centroid Force
    dyn.Nc_ir = Nc_ir; % Centroid Torque
    dyn.f_ir = f_ir; % Output force
    dyn.tau_ir = tau_ir; % Output torque
    dyn.motor_t_ir = motor_t_ir; % Motor output force/torque

    dyn.Fc_l = Fc_l; % Centroid Force
    dyn.f_l = f_l; % Motor output torque
    dyn.tau_l = tau_l; % Motor output force
    dyn.motor_t_l = motor_t_l; % Final motor force/torque
    disp('Dynamics: time ' + string(toc) + 'sec');
    %% Inertial parameter linearization ->tau=H*P
    % Method 1
    tic;
    dyn.H_l = equationsToMatrix(motor_t_l,rbt_df.link_params);
    disp('Linearization Method 1: time ' + string(toc) + 'sec');
    % Method 2
    tic;
    % Auxiliary matrices B A 
    dof = rbt_df.dof;
    B = cell(1,dof);
    A = cell(1,dof);
    for i = 1:dof
        B{i} = getB(w{i+1},dw{i+1},dv{i+1});
        A{i} = getA(w{i+1},dw{i+1},dv{i+1});
    end

    % Auxiliary matrices Yf, Yn
    Yf = cell(1,dof);
    Yn = cell(1,dof);
    Yf{dof}=sym([zeros(3,(dof-1)*10) B{dof}]); Yn{dof}=sym([zeros(3,(dof-1)*10) A{dof}]); 
    for i = dof-1:-1:1
        [Yf{i}, Yn{i}]=get_Yf_Yn(i,geom.Rn_n1{i+1},B{i},A{i},Yf{i+1},Yn{i+1},geom.Pn_n1{i+1},dof);
    end
    % Linear matrix H
    H = sym(zeros(dof,dof*10));
    for i = 1:dof
        H(i,:) = [0,0,1]*Yn{i};
        disp('Generate H(n,:):'+ string(i))
    end

    dyn.H = H;

    disp('Linearization Method 2: time ' + string(toc) + 'sec');
    
    
end

% Linearization function
function K=getK(vec3)
K=[vec3(1)  vec3(2)  vec3(3)    0        0        0;
     0      vec3(1)    0      vec3(2)  vec3(3)    0;
     0        0      vec3(1)    0      vec3(2)  vec3(3)];
end

% Cross product matrix
function S=getS(vec3)
S=[   0     -vec3(3)  vec3(2);
   vec3(3)     0     -vec3(1);
   -vec3(2) vec3(1)     0    ];
end

% Torque linearization matrix H 3x10
function B=getB(w,dw,dv)
B = sym(zeros(3,10));
B(:,7:9)=getS(dw)+getS(w)*getS(w);
B(:,10)=dv;
end

% Force linearization matrix A 3x10
function A=getA(w,dw,dv)
A = sym(zeros(3,10));
A(:,1:6)=getK(dw)+getS(w)*getK(w);
A(:,7:9)=-getS(dv);
end

function [Yf, Yn] = get_Yf_Yn(id, R, H, A, Yf_next, Yn_next, Po,dof)
c1=(id-1)*10+1;
c2=(id-1)*10+10;

Yf=sym(zeros(3,dof*10));
Yf(:,c1:c2)=H;
Yf=Yf+R*Yf_next;

Yn=sym(zeros(3,dof*10));
Yn(:,c1:c2)=A;
Yn=Yn+R*Yn_next+getS(Po)*R*Yf_next;
end