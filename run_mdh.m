clear
clc
%%
gravity_vec = [0; 0; 9.80665];

rbt_config = Demo_MDH;
rbt = struct();
rbt.rbt_df = DefineRobot_MDH('Demo_MDH',rbt_config);
rbt.geom = GeometryCalculation_MDH(rbt.rbt_df);
rbt.dyn = Dynamics_MDH(rbt.rbt_df, rbt.geom, gravity_vec);
rbt.base = DynamicBaseParamCalc_MDH(rbt.rbt_df, rbt.dyn);


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