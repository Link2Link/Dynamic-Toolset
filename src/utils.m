function result = utils
    result.vec2so3 = @vec2so3;
    result.so32vec = @so32vec;
    result.vec2se3 = @vec2se3;
    result.se32vec = @se32vec;
    result.matrix_exp_3 = @matrix_exp_3;
    result.matrix_exp_6 = @matrix_exp_6;
    result.inertia_vec2tensor = @inertia_vec2tensor;
    result.inertia_tensor2vec = @inertia_tensor2vec;
    result.translation_transfmat = @translation_transfmat;
    result.ml2r = @ml2r;
    result.Lmr2I = @Lmr2I;
    result.gen_DLki_mat = @gen_DLki_mat;
    result.gen_DLki_mat4 = @gen_DLki_mat4;
    result.PhyConsMat = @PhyConsMat;
end

function matrix = vec2so3(vector)
    matrix = [0, -vector(3), vector(2); 
        vector(3), 0, -vector(1);
        -vector(2), vector(1), 0];
end

function vector = so32vec(matrix)
    vector = [matrix(3, 2);
        matrix(1, 3);
        matrix(2, 1)];
end

function matrix = vec2se3(vector)
    matrix = [vec2so3(vector(1:3)), vector(4:6);
        0, 0, 0, 0];
end

function vector = se32vec(matrix)
    vector = [matrix(3, 2); 
        matrix(1, 3); 
        matrix(2, 1); 
        matrix(1:3, 4)];
end

function  R = matrix_exp_3(omghat, theta)
    omg_mat = vec2so3(omghat);
    R = eye(3) + sin(theta)*omg_mat + (1-cos(theta))*omg_mat*omg_mat;
end

function T = matrix_exp_6(screw, theta)
    omghat = screw(1:3);
    screw_mat = vec2se3(screw);
    
    omg_mat = screw_mat(1:3, 1:3);
    T = [matrix_exp_3(omghat, theta), (eye(3)*theta + (1-cos(theta))*omg_mat + (theta-sin(theta))*omg_mat*omg_mat)*screw_mat(1:3, 4);
         0, 0, 0, 1];
end

function I = inertia_vec2tensor(vector)
    I = [vector(1), vector(2), vector(3);
        vector(2), vector(4), vector(5);
        vector(3), vector(5), vector(6)];
end

function vector = inertia_tensor2vec(I)
    vector = [I(1, 1), I(1, 2), I(1, 3), I(2, 2), I(2, 3), I(3, 3)];
end

function matrix = translation_transfmat(v)
    matrix = [1, 0, 0, v(1);
        0, 1, 0, v(2);
        0, 0, 1, v(3);
        0, 0, 0, 1];
end

function r = ml2r(m, l)
    r = l / m;
end

function I = Lmr2I(L, m, r)
    I = L - m * vec2so3(r)' * vec2so3(r);
end

function M = gen_DLki_mat()
    M = cell(1, 10);
    for i = 1:10
        M{i} = cell(6, 6);
    end
    % Lxx
    M{1}{1, 1} = 1;
    % Lxy
    M{2}{1, 2} = 1;
    M{2}{2, 1} = 1;
    % Lxz
    M{3}{1, 3} = 1;
    M{3}{3, 1} = 1;
    % Lyy
    M{4}{2, 2} = 1;
    % Lyz
    M{5}{2, 3} = 1;
    M{5}{3, 2} = 1;
    % Lzz
    M{6}{3, 3} = 1;
    % lx
    M{7}{2, 6} = 1;
    M{7}{6, 2} = 1;
    M{7}{3, 5} = -1;
    M{7}{5, 3} = -1;
    % ly
    M{8}{1, 6} = -1;
    M{8}{6, 1} = -1;
    M{8}{3, 4} = 1;
    M{8}{4, 3} = 1;
    % lz
    M{9}{1, 5} = 1;
    M{9}{5, 1} = 1;
    M{9}{2, 4} = -1;
    M{9}{4, 2} = -1;
    % m
    M{10}{4, 4} = 1;
    M{10}{5, 5} = 1;
    M{10}{6, 6} = 1;
end

function M = gen_DLki_mat4()
    M = cell(1, 10);
    for i = 1:10
        M{i} = cell(4, 4);
    end
    % Lxx
    M{1}{1, 1} = -0.5;
    M{1}{2, 2} = 0.5;
    M{1}{3, 3} = 0.5;

    % Lxy
    M{2}{1, 2} = -1;
    M{2}{2, 1} = -1;
    % Lxz
    M{3}{1, 3} = -1;
    M{3}{3, 1} = -1;
    % Lyy
    M{4}{1, 1} = 0.5;
    M{4}{2, 2} = -0.5;
    M{4}{3, 3} = 0.5;

    % Lyz
    M{5}{2, 3} = -1;
    M{5}{3, 2} = -1;
    % Lzz
    M{6}{1, 1} = 0.5;
    M{6}{2, 2} = 0.5;
    M{6}{3, 3} = -0.5;
    % lx
    M{7}{1, 4} = 1;
    M{7}{4, 1} = 1;
    % ly
    M{8}{2, 4} = 1;
    M{8}{4, 2} = 1;
    % lz
    M{9}{3, 4} = 1;
    M{9}{4, 3} = 1;
    % m
    M{10}{4, 4} = 1;
end

function D = PhyConsMat(delta)
    %delta = [Lxx, Lxy, Lxz, Lyy, Lyz, Lzz, lx, ly, lz, m]
    L = inertia_vec2tensor(delta(1:6));
    l = delta(7:9);
    m = delta(10);
    D0 = trace(L)/2*eye(3)-L;
    D = [D0, l; l', m];
end