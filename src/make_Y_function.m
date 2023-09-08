function func_handle = make_Y_function(rbt)

out_path = "output/";
[status, msg, msgID] = mkdir('output');

out_path = "output/";
Y = rbt.base.H_b;
func_handle = matlabFunction(Y, 'File',out_path+rbt.rbt_df.name+'_Y_func', 'Vars', {rbt.rbt_df.coordinates', rbt.rbt_df.d_coordinates', rbt.rbt_df.dd_coordinates'});
disp('function writen into '+ out_path+rbt.rbt_df.name + '_Y_func.m')

matlabFunction(vpa(rbt.base.prime2beta), 'File',out_path+rbt.rbt_df.name+'_prime_param2beta_param');
disp('prime param to beta param mapping matrix writen into '+ out_path + rbt.rbt_df.name + '_prime_param2beta_param.m')


end