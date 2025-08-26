function rbt_config = DemoTianyiArmR()

%% 7 joint angle
syms q1 q2 q3 q4 q5 q6 q7 real;
pi_ = sym(pi);

%% geoetric parameter
R = eye(3);
p{1} = [0 0 0]';                %零位状态下左手其他坐标系到根部的位移
p{2} = [0 -0.068      0]';
p{3} = [0  -0.068 -0.1025]';
p{4} = [0.02 -0.068   -0.3]';
p{5} = [0  -0.068 -0.3543]';
p{6} = [0  -0.07 -0.543]';
p{7} = [0  -0.07 -0.563]';

ws = [0 1 0
1 0 0
0 0 1
0 1 0
0 0 1
0 1 0
1 0 0]';

% define link number used to create joint tree
L_b = 0;
L_1 = 1;
L_2 = 2;
L_3 = 3;
L_4 = 4;
L_5 = 5;
L_6 = 6;
L_7 = 7;

% define spring
dlN = NaN;

rbt_config = {
% Joint number|prev link|succ links     |screw axis unit direction  |rs                 |M_R                    |theta      |h    |link inertia     |spring 
    L_b,          -1,       L_1,         [0;0;0],                  [0;0;0],             eye(3),                sym(0),     0,     false,         dlN;
    L_1,          L_b,      L_2,         ws(:, 1),                 p{1},                eye(3),                q1,         0,     true,          dlN;
    L_2,          L_1,      L_3,         ws(:, 2),                 p{2},                eye(3),                q2,         0,     true,          dlN; 
    L_3,          L_2,      L_4,         ws(:, 3),                 p{3},                eye(3),                q3,         0,     true,          dlN;
    L_4,          L_3,      L_5,         ws(:, 4),                 p{4},                eye(3),                q4,         0,     true,          dlN;
    % L_5,          L_4,      L_6,         ws(:, 5),                 p{5},                eye(3),                q5,         0,     true,          dlN; 
    % L_6,          L_5,      L_7,         ws(:, 6),                 p{6},                eye(3),                q6,         0,     true,          dlN;
    % L_7,          L_6,      [],          ws(:, 7),                 p{7},                eye(3),                q7,         0,     true,          dlN
};
end
