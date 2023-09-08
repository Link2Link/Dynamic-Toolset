function result = make_MCG_function(rbt)
out_path = "output/";
[status, msg, msgID] = mkdir('output');


tic
disp('Creating M C G matrix ...');
model_name = rbt.rbt_df.name;
P = sym('p', [rbt.base.base_num,1], 'real');
reresult.base_param_symbol = P;

DOF = size(rbt.base.H_b, 1);
Y = rbt.base.H_b;
tau = Y*P;

%% seperate m g c
g = subs(tau, [rbt.rbt_df.d_coordinates, rbt.rbt_df.dd_coordinates], [zeros(size(rbt.rbt_df.d_coordinates)), zeros(size(rbt.rbt_df.dd_coordinates))]);
mc = tau - g;
m = subs(mc, [rbt.rbt_df.d_coordinates], [zeros(size(rbt.rbt_df.d_coordinates))]);
c = mc - m;
reresult.g = g;
reresult.m = m;
reresult.c = c;


%% make static part
vars = symvar(g);
static_part_vars = [];
for v = vars
    if ~ismember(v, rbt.rbt_df.coordinates)
        static_part_vars = [static_part_vars; v];
    end
end
[base2static_part, ~] = equationsToMatrix(static_part_vars, P);
reresult.static_param_symbol = static_part_vars;
reresult.base_param2static_param = base2static_part;
[Gmat, ~] = equationsToMatrix(g, static_part_vars);
reresult.Gmat = simplify(Gmat);

matlabFunction(reresult.Gmat, 'File',out_path+model_name+'_G_mat', 'Vars', {rbt.rbt_df.coordinates'});
disp('function writen into '+ out_path + model_name + '_G_mat.m')

matlabFunction(base2static_part, 'File',out_path+model_name+'_base_param2static_param');
disp('base param to static param mapping matrix writen into '+ out_path + model_name + '_base_param2static_param.m')

%% make G matrix
G = vpa(simplify(g));
matlabFunction(G, 'File',out_path+model_name+'_G_func', 'Vars', {rbt.rbt_df.coordinates', P});
disp('function writen into '+ out_path + model_name + '_G_func.m')


%% make M matrix
M = [];
for k = 1:DOF
    one_hot = zeros(DOF, 1);
    one_hot(k) = 1;
    M_col = subs(m, [rbt.rbt_df.dd_coordinates], [one_hot']);
    M = [M, M_col];
end
matlabFunction(vpa(M), 'File',out_path+model_name+'_M_func', 'Vars', {rbt.rbt_df.coordinates', P});
disp('function writen into '+ out_path + model_name + '_M_func.m')


%% make C matrix
q = rbt.rbt_df.coordinates;
dq = rbt.rbt_df.d_coordinates;
Christoffel = @(i,j,k) 0.5*(diff(M(i,j), q(k)) + diff(M(i,k), q(j)) - diff(M(j,k), q(i)));
C = sym(zeros(DOF));
for i = 1:DOF
    for j = 1:DOF
        C(i,j) = 0;
        for k = 1:DOF
            C(i,j) = C(i,j) + Christoffel(i,j,k)*dq(k);
        end
    end
end
rbt.C_b_func = matlabFunction(vpa(C), 'File',out_path+model_name+'_C_func', 'Vars', {rbt.rbt_df.coordinates', rbt.rbt_df.d_coordinates', P});
disp('function writen into '+ out_path + model_name + '_C_func.m')

%% check
ddq = rbt.rbt_df.dd_coordinates;
Tau = M*ddq(:) + C*dq(:) + g;
err = Tau - tau;
rand_err = subs(err, [q, dq, ddq, P'], [rand(size(q)), rand(size(dq)), rand(size(ddq)), rand(size(P'))]);
abs_err = vpa(norm(rand_err));
disp("Symbol deduce error is : "+string(abs_err))

time_consum = toc;
disp("make_MCG_func takes " + string(time_consum) + " sec");


end