clear
clc
%%
gravity_vec = [0, 0, -9.8015];

rbt_config = Demo7DOF;
rbt = struct();
rbt.rbt_df = DefineRobot('Demo7DOF',rbt_config);
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