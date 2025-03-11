function [traj, result] = trajactory_optimization(rbt, qmax, qmin, dqmax, ddqmax, L, Tf, N_pop, N_iter)
dim = length(qmax);
result = struct();
result.cost = [];
%% 参数存储
traj = struct();
traj.dim = dim;
traj.L = L;                        %傅里叶级数
traj.Tf = Tf;                      %轨迹周期 sec
traj.wf = 2*pi/Tf;
traj.qmax = qmax(:);
traj.qmin = qmin(:);
traj.dqmax = dqmax(:);
traj.ddqmax = ddqmax(:);

Y_handle = matlabFunction(rbt.base.H_b, 'Vars', {rbt.rbt_df.coordinates', rbt.rbt_df.d_coordinates', rbt.rbt_df.dd_coordinates'});

%% 生成轨迹符号表达式
a = sym("a", [traj.dim, traj.L], "real");
b = sym("b", [traj.dim, traj.L], "real");
q0_b = sym("q0_b", [traj.dim, 1], "real");
dq0_b = sym("dq0_b", [traj.dim, 1], "real");
ddq0_b = sym("ddq0_b", [traj.dim, 1], "real");
t = sym("t", "real");
wf = traj.wf;
dt = sym("dt", "real");
% qbias = sym("qb", [traj.dim, 1]);


% calc base trajactory
for k = 1:traj.dim
    phy_base{k} = sym(0);
    for i = 1:traj.L
        phy_base{k} = phy_base{k} + a(k, i) / i / wf * sin(i*wf*t);
        phy_base{k} = phy_base{k} - b(k, i) / i / wf * cos(i*wf*t);
    end
end

% calc offset trajactory
for k = 1:traj.dim
    phy_offset{k} = q0_b(k) - dq0_b(k)/wf*sin(wf*t) + ddq0_b(k)/wf/wf*cos(wf*t);
end

% full trajactory
for k = 1:traj.dim
    traj.phy(k, 1) = phy_base{k} + phy_offset{k};
    dphy_dt = TimeDerivative(traj.phy(k, 1), t, 1);
    ddphy_ddt = TimeDerivative(dphy_dt, t, 1);
    traj.speed(k, 1) = dphy_dt;
    traj.acc(k, 1) = ddphy_ddt;
end

traj.func = matlabFunction(traj.phy, 'Vars', {a, b, q0_b, dq0_b, ddq0_b, t});
traj.func_speed = matlabFunction(traj.speed, 'Vars', {a, b, q0_b, dq0_b, ddq0_b, t});
traj.func_acc = matlabFunction(traj.acc, 'Vars', {a, b, q0_b, dq0_b, ddq0_b, t});
save("output/traj.mat", "traj");
%% TLBO

% 随机初始化
for k = 1:N_pop
    param{k} = struct();
    param{k}.a = (rand(dim, L) - 0.5) .* (traj.qmax - traj.qmin);
    param{k}.b = (rand(dim, L) - 0.5) .* (traj.qmax - traj.qmin);
    param{k}.q0_b = rand(dim, 1);
    param{k}.dq0_b = rand(dim, 1);
    param{k}.ddq0_b = rand(dim, 1);
    param{k} = process_constraint(traj, param{k});
    cost(k) = calc_cost(traj, param{k}, Y_handle);
end

for loop = 1:N_iter
tic;

%% 选取教师
[min_cost, idx] = min(cost);
X_tch = param{idx};
J_tch = cost(k);
result.param = X_tch;
result.cost = [result.cost, J_tch];

save("output/traj_optimizaiton_result.mat", "result");

% 班级平均值
Xmean.a = zeros(dim, L);
Xmean.b = zeros(dim, L);
Xmean.q0_b = zeros(dim, 1);
Xmean.dq0_b = zeros(dim, 1);
Xmean.ddq0_b = zeros(dim, 1);
for k = 1:N_pop
    Xmean.a         = Xmean.a       + param{k}.a / N_pop;
    Xmean.b         = Xmean.b       + param{k}.b / N_pop;
    Xmean.q0_b      = Xmean.q0_b    + param{k}.q0_b / N_pop;
    Xmean.dq0_b     = Xmean.dq0_b   + param{k}.dq0_b / N_pop;
    Xmean.ddq0_b    = Xmean.ddq0_b  + param{k}.ddq0_b / N_pop;
end
Jmean = calc_cost(traj, Xmean, Y_handle);

% 教学策略 Dmean = (X_tch - TF*Xmean) .* ri

TF =  round(1 + rand());
Dmean.a         =     (X_tch.a       -   TF*Xmean.a).*rand(traj.dim, traj.L);
Dmean.b         =     (X_tch.b       -   TF*Xmean.b).*rand(traj.dim, traj.L);
Dmean.q0_b      =     (X_tch.q0_b    -   TF*Xmean.q0_b).*rand(traj.dim, 1);
Dmean.dq0_b     =     (X_tch.dq0_b   -   TF*Xmean.dq0_b).*rand(traj.dim, 1);
Dmean.ddq0_b    =     (X_tch.ddq0_b  -   TF*Xmean.ddq0_b).*rand(traj.dim, 1);

% 学生更新
for k = 1:N_pop
    param_new{k}.a         = param{k}.a       + Dmean.a      ;
    param_new{k}.b         = param{k}.b       + Dmean.b      ;
    param_new{k}.q0_b      = param{k}.q0_b    + Dmean.q0_b   ;
    param_new{k}.dq0_b     = param{k}.dq0_b   + Dmean.dq0_b  ;
    param_new{k}.ddq0_b    = param{k}.ddq0_b  + Dmean.ddq0_b ;
    param_new{k} = process_constraint(traj, param_new{k});  % 新策略也要满足约束
    cost_new = calc_cost(traj, param_new{k}, Y_handle);
    if cost_new < cost(k)   %若新策略更优，则选择新策略
        param{k} = param_new{k};
        cost(k) = cost_new;
    end
end

% 学生相互学习
for d = 1:N_pop
    % 学生d向s学习
    for s = 1:N_pop
        if d == s
            continue
        end
    
        if cost(d) > cost(s)
            X1 = param{s};
            X2 = param{d};
        else
            X1 = param{d};
            X2 = param{s};
        end
        
        % Xd_new = Xd + rd .* (X1 - X2) 
        X_new.a           =   param{d}.a          + (X1.a      -       X2.a        ) .* rand(traj.dim, traj.L);
        X_new.b           =   param{d}.b          + (X1.b      -       X2.b        ) .* rand(traj.dim, traj.L);
        X_new.q0_b        =   param{d}.q0_b       + (X1.q0_b   -       X2.q0_b     ) .* rand(traj.dim, 1);
        X_new.dq0_b       =   param{d}.dq0_b      + (X1.dq0_b  -       X2.dq0_b    ) .* rand(traj.dim, 1);
        X_new.ddq0_b      =   param{d}.ddq0_b     + (X1.ddq0_b -       X2.ddq0_b   ) .* rand(traj.dim, 1);
        X_new = process_constraint(traj, X_new);  % 新策略也要满足约束
        cost_new = calc_cost(traj, X_new, Y_handle);
        if cost_new < cost(k)   %若新策略更优，则选择新策略
            param{k} = X_new;
            cost(k) = cost_new;
        end
    end
end
looptime = toc;
disp(" loop " + string(loop) + " used time " + string(looptime) + " cost = " + string(min(cost)));
end

[min_cost, idx] = min(cost);
X_tch = param{idx};
J_tch = cost(k);
disp('optimization end iter = '+ string(loop) + ', cost = ' +string(J_tch))

end




function cost = calc_cost(traj, param, Y_handle)
W = [];
% ts = 0:(2*traj.L + 1):traj.Tf;

for t =  0:(traj.Tf/(2*traj.L + 1)):traj.Tf
    q = traj.func(param.a, param.b, param.q0_b, param.dq0_b, param.ddq0_b, t);
    dq = traj.func_speed(param.a, param.b, param.q0_b, param.dq0_b, param.ddq0_b, t);
    ddq = traj.func_acc(param.a, param.b, param.q0_b, param.dq0_b, param.ddq0_b, t);
    Y = Y_handle(q, dq, ddq);
    W = [W; Y];
end
cost = cond(W);
end



%% 约束处理
function param = process_constraint(traj, param)
    dim = traj.dim;
    L = traj.L;
    a = param.a;
    b = param.b;
    q0_b = param.q0_b;
    dq0_b = param.dq0_b;
    ddq0_b = param.ddq0_b;
    wf = traj.wf;

    for k = 1:dim 
        dq0_b(k) = 0;
        ddq0_b(k) = 0;
        for i = 1:L
            dq0_b(k) = dq0_b(k) + a(k, i);
            ddq0_b(k) = ddq0_b(k) + i * wf * b(k, i);
        end
    end

    num = 10*L;
    ts = 0:traj.Tf/num:traj.Tf;
    q_sam = zeros(dim, num);
    dq_sam = zeros(dim, num);
    ddq_sam = zeros(dim, num);
    for k = 1 : num
        q_sam(:, k) = traj.func(a, b, q0_b, dq0_b, ddq0_b, ts(k));
        dq_sam(:, k) = traj.func_speed(a, b, q0_b, dq0_b, ddq0_b, ts(k));
        ddq_sam(:, k) = traj.func_acc(a, b, q0_b, dq0_b, ddq0_b, ts(k));
    end
    QMax = max(q_sam, [], 2);
    QMin = min(q_sam, [], 2);
    dQ = max(abs(dq_sam), [], 2);
    ddQ = max(abs(ddq_sam), [], 2);
    
    lambda_v = traj.dqmax./dQ;                              % 速度缩放因子
    lambda_a = traj.ddqmax./ddQ;                            % 加速度缩放因子
    lambda_p = (traj.qmax - traj.qmin) ./ (QMax - QMin);    % 位置缩放因子
    
    lambda = min([lambda_p, lambda_v, lambda_a], [], 2);

    a = lambda .* a;
    b = lambda .* b;
    q0_b = lambda .* q0_b;
    dq0_b = lambda .* dq0_b;
    ddq0_b = lambda .* ddq0_b;

    q_cc = (traj.qmax + traj.qmin) / 2;
    Q_cc = (lambda .* QMax + lambda .* QMin) / 2;
    delta = Q_cc - q_cc;
    q0_b = q0_b - delta;

    param.a = a;
    param.b = b;
    param.q0_b = q0_b;
    param.dq0_b = dq0_b;
    param.ddq0_b = ddq0_b;
end

function df = TimeDerivative(f, q, varargin)
    order = length(varargin);
    
    
    syms t real;
    for k = 1:length(q)
        qt(k) = str2sym(string(q(k))+"t(t)");
        assume(qt(k), 'real');
    end
    
    for k = 1:order
        dq{k} = varargin{k};
        dqt{k} = diff(qt, t, k);
    end
    
    ft = subs(f, q, qt);
    dft = diff(ft, t, order);
    df = dft;
    for k = order:-1:1
        df = subs(df, dqt{k}, dq{k});
    end
    
    
    df = subs(df, qt, q);


end
