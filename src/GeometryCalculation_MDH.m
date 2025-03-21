function geom = GeometryCalculation_MDH(rbt_df)
    % Generate kinematics of robot
    tic;
    format long;
    Tn_n1 = cell(1,rbt_df.dof);
    Rn_n1 = cell(1,rbt_df.dof);
    Pn_n1 = cell(1,rbt_df.dof);

    for i=1:rbt_df.dof
        alpha = rbt_df.alpha{i+1}; 
        a = rbt_df.a{i+1}; 
        d = rbt_df.d{i+1}; 
        theta = rbt_df.theta{i+1} + rbt_df.offset{i+1}; 
        
        Tn_n1{i} = mdh_transform(alpha, a, d, theta);
        Rn_n1{i} = Tn_n1{i}(1:3,1:3);  % Extract the rotation matrix
        Pn_n1{i} = Tn_n1{i}(1:3,4);    % Extract the displacement matrix
    end

    geom.Tn_n1 = Tn_n1;
    geom.Rn_n1 = Rn_n1;
    geom.Pn_n1 = Pn_n1;

    disp('GeometryCalculation: time ' + string(toc) + 'sec');
end

function T = mdh_transform(alpha, a, d, theta)
    if alpha==0
        Salpha=0;
        Calpha=1;
    elseif abs(alpha)==pi/2
        Calpha=0;
        if alpha==pi/2
            Salpha=1;
        elseif alpha==-pi/2
            Salpha=-1;
        end
    elseif abs(alpha)==pi
        Salpha=0;
        Calpha=-1;
    else
        Salpha=sin(alpha);
        Calpha=cos(alpha);
    end
    T = [cos(theta)         -sin(theta)         0          a;
         sin(theta)*Calpha  cos(theta)*Calpha   -Salpha    -d*Salpha;
         sin(theta)*Salpha  cos(theta)*Salpha   Calpha     d*Calpha;
         0                  0                   0          1];
    % Modified DH Transform Matrix (Order: Rot_x -> Trans_x -> Rot_z -> Trans_z)
    % T = [cos(theta)        -sin(theta)          0         a;
    %      sin(theta)*cos(alpha) cos(theta)*cos(alpha) -sin(alpha) -sin(alpha)*d;
    %      sin(theta)*sin(alpha) cos(theta)*sin(alpha)  cos(alpha)  cos(alpha)*d;
    %      0                   0                  0         1];
end