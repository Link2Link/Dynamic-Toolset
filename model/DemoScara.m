function rbt_config = DemoScara()

%% 7 joint angle
syms q1 q2 q3 q4 real;
pi_ = sym(pi);

%% geoetric parameter
d1 = 0.1;
d2 = 0.1;
d3 = 0.1;

% define link number used to create joint tree
L_b = 0;
L_1 = 1;
L_2 = 2;
L_3 = 3;
L_4 = 4;

% define spring
dlN = NaN;

dl4 = [sin(q4), cos(q4), 1];

rbt_config = {
% Joint number|prev link|succ links     |screw axis unit direction  |rs                 |M_R                    |theta      |h    |link inertia     |spring 
    L_b,          -1,       L_1,         [0;0;0],                  [0;0;0],             eye(3),                sym(0),       0,     false,         dlN;
    L_1,          L_b,      L_2,         [0;0;1],                  [0;0;0],             eye(3),                q1,           0,     true,          dlN;
    L_2,          L_1,      L_3,         [0;0;1],                  [d1;0;0],            eye(3),                q2,           0,     true,          dlN; 
    L_3,          L_2,      L_4,         [0;0;1],                  [d1+d2;0;0],         eye(3),                q3,           0,     true,          dlN;
    L_4,          L_3,      [],         [0;0;1],                   [d1+d2;0;-d3],       eye(3),                q4,           inf,   true,          dl4;
};
end
