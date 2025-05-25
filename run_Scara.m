clear
clc
close all 
%%
gravity_vec = [0, 0, -9.8015];

rbt_config = DemoScara;
rbt = struct();
rbt.rbt_df = DefineRobot('DemoScara',rbt_config);
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
result_MCG = make_MCG_function(rbt); 

%% 辨识轨迹优化

qmax = deg2rad([90, 90, 90, 90, ]);
qmin = deg2rad([-90, -90, -90,-90]);
dqmax = [10, 10, 10, 10];
ddqmax = dqmax * 3;

L = 7;              % 傅里叶级数
Tf = 30;            % 周期
N_pop = 100;        % 种群数
N_iter = 20;       % 迭代次数

profile clear;
profile off;
profile on;
[traj, result] = trajactory_optimization(rbt, qmax, qmin, dqmax, ddqmax, L, Tf, N_pop, N_iter);
profile viewer;
profile off;
%% draw traj 
draw_trajactory(traj, result.param)
figure
semilogy(result.cost);
grid on;
xlabel('iter');
ylabel('log(cost)');
