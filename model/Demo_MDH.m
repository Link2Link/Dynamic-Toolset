function robot_cfg = Demo_MDH()

%% geoetric parameter
% define link number used to create joint tree
L_b = 0;
L_1 = 1;
L_2 = 2;
L_3 = 3;
L_4 = 4;
L_5 = 5;
L_6 = 6;

% offset
L0 = 0.07;
L1 = 0.14;
L2 = 0.12;
L3 = 0.051;
L4 = 0.051;

%% Kinematics Parameters
syms q1 q2 q3 q4 q5 q6 real;
mdh = {
% a       |  alpha   |    d      |theta    |offset
0,          0,          L0,         q1,        0;
0,          pi/2,        0,         q2,     pi/2;
L1,         pi,          0,         q3,        0;
L2,        -pi,         L3,         q4,    -pi/2;
0,         -pi/2,       L4,         q5,        0;
0,          pi/2,        0,         q6,        0;
};

% Link Type 0=Revolute 1=Assitive 2=Prismatic 3=Fixed

robot_cfg = {
% Joint number|prev link|succ links     |a         |alpha      |d        |theta      |offset    |Link Type  |use interia
    L_b,          -1 ,      L_1,         0,             0,          0,          0,          0,      3,       false;
    L_1,          L_b,      L_2,         mdh{1,1},   mdh{1,2},  mdh{1,3},  mdh{1,4},   mdh{1,5},    0,       true;
    L_2,          L_1,      L_3,         mdh{2,1},   mdh{2,2},  mdh{2,3},  mdh{2,4},   mdh{2,5},    0,       true;
    L_3,          L_2,      L_4,         mdh{3,1},   mdh{3,2},  mdh{3,3},  mdh{3,4},   mdh{3,5},    0,       true;
    L_4,          L_3,      L_5,         mdh{4,1},   mdh{4,2},  mdh{4,3},  mdh{4,4},   mdh{4,5},    0,       true;
    L_5,          L_4,      L_6,         mdh{5,1},   mdh{5,2},  mdh{5,3},  mdh{5,4},   mdh{5,5},    0,       true; 
    L_6,          L_5,      [],          mdh{6,1},   mdh{6,2},  mdh{6,3},  mdh{6,4},   mdh{6,5},    0,       true;
};