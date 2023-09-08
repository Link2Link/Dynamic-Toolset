function func_handle = make_H_function(rbt)

out_path = "output/";
[status, msg, msgID] = mkdir('output');

out_path = "output/";
H = rbt.base.H;
func_handle = matlabFunction(H, 'File',out_path+rbt.rbt_df.name+'_H_func', 'Vars', {rbt.rbt_df.coordinates', rbt.rbt_df.d_coordinates', rbt.rbt_df.dd_coordinates'});
disp('function writen into '+ out_path+rbt.rbt_df.name + '_H_func.m')

end