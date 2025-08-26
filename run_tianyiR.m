clear
clc
%%
syms gx gy gz real;
gravity_vec = [0, 0, -9.8015];
gravity_vec = [gx, gy, gz];

rbt_config = DemoTianyiArmR;
rbt = struct();
rbt.rbt_df = DefineRobot('DemoTianyiArmR',rbt_config);
rbt.geom = GeometryCalculation(rbt.rbt_df);
rbt.dyn = Dynamics(rbt.rbt_df, rbt.geom, gravity_vec);
rbt.base = DynamicBaseParamCalc(rbt.rbt_df, rbt.dyn);


%% save output
out_path = "output/";
[status, msg, msgID] = mkdir('output');
addpath(out_path);

save output/rbt;

%% H matrix for prime parameter
make_H_function(rbt);

%% Y matrix for base parameter 
make_Y_function(rbt);

%% MCG matrix : M(q, P)  C(q, dq, P)  G(q, P) 
make_MCG_function(rbt); 

% %% trajactory optimization
% qmax = deg2rad([90, 90, 90, 90, 90, 90, 180]);
% qmin = deg2rad([-90, -90, -90, -90, -90, -90, -180]);
% dqmax = [10, 10, 10, 10, 10, 10, 10];
% ddqmax = dqmax * 3;
% 
% L = 7;              % 傅里叶级数
% Tf = 30;            % 周期
% N_pop = 100;        % 种群数
% N_iter = 20;       % 迭代次数
% 
% [traj, result] = trajactory_optimization(rbt, qmax, qmin, dqmax, ddqmax, L, Tf, N_pop, N_iter);
% 
% %% draw traj 
% draw_trajactory(traj, result.param)
% figure
% semilogy(result.cost);
% grid on;
% xlabel('iter');
% ylabel('log(cost)');
