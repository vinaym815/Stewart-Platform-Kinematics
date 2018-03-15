clear all;
close all;
clc;

%% 
% Each column of legs have the structure [px; py; pz; bx, by, bz]; where p
% stands for platform and b stands for base. (px, py, pz) are the coordinates
% of the point on leg attached to platform, represented in the frame 
% rigidly attached to the center of top platform. Similarly (bx, by, bz) are
% the coordinated of point on leg attached to base, represented in the frame
% rigidly attached to the base centroid.

legs = [100,    50,     -87,    -100,   -50,    87;
        0,      87,     50,     0,      -87,    -50;
        0,      0,      0,      0,      0,      0;
        150,    75,     -130,   -150,   -75,    130;
        0,      130,    75,     0,      -130,   -75;
        0,      0,      0,      0,      0,      0];

% Pf = [Yaw, Pitch, Roll, X, Y, Z]  platform information represented in fixed frame 
Pf = [degtorad(0); degtorad(0); degtorad(0); 0; 0; 300];

%%
% Transformation matrix representing platform orientation wrt base
% Following the convention yaw pitch roll
thz = Pf(1);    thy = Pf(3);    thx = Pf(2);
Rz = [cos(thz)      -sin(thz)   0;
     sin(thz)       cos(thz)    0;
     0              0           1];

Ry = [cos(thy)      0   sin(thy);
     0              1          0;
     -sin(thy)      0   cos(thy)];

Rx = [1     0           0;
     0      cos(thx)    -sin(thx);
     0      sin(thx)    cos(thx)];

Rot = Rz*Ry*Rx;

PlatTran = [Rot, [Pf(4); Pf(5); Pf(6)];
            zeros(1,3), 1];
%%
% Initializing leg lengths
LegLengths = zeros(6,1);
jac = zeros(6,6);   % Kinematic Jacobian about the origin of platform
% As the actuators are applying only forces P*, Q*, R* simplifies to P, Q, R
% (P, Q, R), (P*, Q*, R*) and (L, M, N) notations are from the screw theory

% Iterating over each leg
for i = 1:6
    % Representing the platform attachement point of leg in coordinate of base
    % platform
    PlatPt = PlatTran*[legs(1:3,i);1];          

    % r is the vector from the base of the platform to leg attachment point on
    % the platform represented in base frame
    r = PlatPt(1:3) - Pf(4:6);

    BasePt = legs(4:6,i);
    
    L = PlatPt(1) - BasePt(1);  M = PlatPt(2) - BasePt(2);  N = PlatPt(3) - BasePt(3) ;
    LegLengths(i) = sqrt(L^2+M^2+N^2);

    jac(i,:) = [L,M,N,cross(r,[L;M;N])']/LegLengths(i);

end

disp(LegLengths)
disp(jac)
