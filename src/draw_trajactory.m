function draw_trajactory(traj, param)
    ts = 0:1e-2:traj.Tf;
    pos = zeros(traj.dim, length(ts));
    for k = 1:length(ts)
        pos(:, k) = traj.func(param.a, param.b, param.q0_b, param.dq0_b, param.ddq0_b, ts(k));
    end
    plot(pos')
    grid on 
end 
