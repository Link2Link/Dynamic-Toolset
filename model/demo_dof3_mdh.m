function robot_cfg = demo_dof3_mdh()

%% geoetric parameter
% define link number used to create joint tree
L_b = 0;
L_1 = 1;
L_2 = 2;
L_3 = 3;

%% Kinematics Parameters
% offset
L0 = 0.07;
L1 = 0.14;

% a
a(1) = 0;
a(2) = 0;
a(3) = L1;

% alpha
alpha(1) = 0;
alpha(2) = pi/2;
alpha(3) = pi;

% d
d(1) = L0;
d(2) = 0;
d(3) = 0;

% theta
syms q1 q2 q3 real;



%% Dynamic Parameters
% Link Type 0=Revolute 1=Assitive 2=Prismatic 3=Fixed

robot_cfg = {
% Joint number|prev link|succ links     |a       |alpha       |d      |theta     |offset |Joint Type  |use interia  
    L_b,          -1,       L_1,         0,      0,            0,       sym(0)      0,          3,         false,;
    L_1,          L_b,      L_2,         a(1),   alpha(1),     d(1),     q1         0,          0,         true,;
    L_2,          L_1,      L_3,         a(2),   alpha(2),     d(2),     q2       pi/2,         0,         true,;
    L_3,          L_2,      [],          a(3),   alpha(3),     d(3),     q3         0,          0,         true,;
};