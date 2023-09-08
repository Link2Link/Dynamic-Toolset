function rbt_config = Demo7DOF()

%% 7 joint angle
syms q1 q2 q3 q4 q5 q6 q7 real;
pi_ = sym(pi);

%% geoetric parameter
d1 = 0.241;
d2 = 0.280;
d3 = 0.330;
d3prime = 0.07;
d4 = 0.1595;

% define link number used to create joint tree
L_b = 0;
L_1 = 1;
L_2 = 2;
L_30 = 3;
L_31 = 4;
L_32 = 5;
L_4 = 6;
L_5 = 7;
L_6 = 8;
L_7 = 9;

% define spring
dlN = NaN;
dl2 = [sin(q2), cos(q2), 1];
dl3 = [sin(q3), cos(q3), 1];
dl5 = [sin(q5), cos(q5), 1];


rbt_config = {
% Joint number|prev link|succ links     |screw axis unit direction  |rs                 |M_R                    |theta      |h    |link inertia     |spring 
    L_b,          -1,       L_1,         [0;0;0],                  [0;0;0],                eye(3),                sym(0),     0,     false,         dlN;
    L_1,          L_b,      [L_2, L_31], [0;0;1],                  [0;0;0],                eye(3),                q1,         0,     true,          dlN;
    L_2,          L_1,      L_30,        [0;-1;0],                 [0;0;-d1],              eye(3),                q2,         0,     true,          dl2; 
    L_30,         L_2,      L_4,         [0;-1;0],                 [0;0;-(d1+d2)],         eye(3),                q3-q2,      0,     true,          dlN;
    L_31,         L_1,      L_32,        [0;-1;0],                 [0;0;-d1],              eye(3),                q3,         0,     true,          dl3; 
    L_32,         L_31,     [],          [0;-1;0],                 [d3prime;0;-d1],        eye(3),                q2-q3,      0,     true,          dlN;
    L_4,          L_30,     L_5,         [0;0;1],                  [d3;0;-(d1+d2-d4)],     eye(3),                q4,         0,     true,          dlN;
    L_5,          L_4,      L_6,         [1;0;0],                  [d3;0;-(d1+d2-d4)],     eye(3),                q5,         0,     true,          dl5; 
    L_6,          L_5,      L_7,         [0;0;1],                  [d3;0;-(d1+d2-d4)],     eye(3),                q6,         0,     true,          dlN;
    L_7,          L_6,      [],          [-1;0;0],                 [d3;0;-(d1+d2-d4)],     [-1,0,0;0,-1,0;0,0,1], q7,         0,     true,          dlN
};
end
