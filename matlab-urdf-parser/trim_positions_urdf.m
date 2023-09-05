function model = trim_positions_urdf(model,VERBOSE)
%
% Make sure all joint positions remain within the range
%
for i = 1:model.n_joint
    lower_limit = model.joint_limits{i}(1);
    upper_limit = model.joint_limits{i}(2);
    joint_name = model.joint_names{i};
    curr_pos = model.positions_rad(i);
    if curr_pos < lower_limit
        if VERBOSE
            fprintf(2,'[trim_positions_within_limits] [%s] [%+4.2f] should be bigger than [%+4.2f].\n',...
                joint_name,curr_pos,lower_limit);
        end
        curr_pos = lower_limit;
    end
    if curr_pos > upper_limit
        if VERBOSE
            fprintf(2,'[trim_positions_within_limits] [%s] [%+4.2f] should be smaller than [%+4.2f].\n',...
                joint_name,curr_pos,upper_limit);
        end
        curr_pos = upper_limit;
    end
    model.positions_rad(i) = curr_pos;
end
