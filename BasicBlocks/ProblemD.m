close all;
clear all;
clc;

% This script does the forward kinematic by jacobian approximation and
% error correction. We don't have a control on the amount of error

%%
% InPfOri = [Roll, Pitch, Yaw, X, Y, Z] initial platform
% position

InPfOri = [degtorad(0); degtorad(20); degtorad(20); 0; 0; 300];

% Platform and base end of legs, in plarform and base coordinates
% respectivelty   [leg1, leg2 ..... leg6]

legs = [100,    50,     -87,    -100,   -50,    87;
        0,      87,     50,     0,      -87,    -50;
        0,      0,      0,      0,      0,      0;
        150,    75,     -130,   -150,   -75,    130;
        0,      130,    75,     0,      -130,   -75;
        0,      0,      0,      0,      0,      0];
%%
thz = InPfOri(1);    thy = InPfOri(2);    thx = InPfOri(3);
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

% InPlatPose is platform pose matrix 
InPlatPose = [Rot, [InPfOri(4); InPfOri(5); InPfOri(6)];
            zeros(1,3), 1];
%%
InLegLn = zeros(6,1);
InJac = zeros(6,6);   
% Kinematic Jacobian about the origin of platform
% As the actuators are applying only forces P*, Q*, R* simplifies to P, Q, R  

% PlatPt is coordinate of leg end attached to platform, represented
% in base frame

% r is the vector joining platform origin to PlatPt represented in
% base frame
for i = 1:6
    PlatPt = InPlatPose*[legs(1:3,i);1];          
    BasePt = legs(4:6,i);
    r = PlatPt(1:3) - InPlatPose(1:3,4);
    
    L = PlatPt(1) - BasePt(1);  M = PlatPt(2) - BasePt(2);  N = PlatPt(3) - BasePt(3) ;
    InLegLn(i) = sqrt(L^2+M^2+N^2);
    InJac(i,:) = [L,M,N,cross(r,[L;M;N])']/InLegLn(i);
end

if rank(InJac)<6 ; disp('Intial configuration of platform was singular'); return; end

clearvars -except legs InLegLn InPlatPose InJac 
%%
ThPlatPose = InPlatPose;
ThLegLn = InLegLn;
ThJac = InJac;
step = 0.1;
for i = step:step:(5+step)                   % sought of jacobian control with error correction
    if rem(round((i - step)*(1/step))/(1/step),1) == 0
%         disp(i)
        disp(ThPlatPose)
        disp(ThLegLn)
    end
    multiplier = 1+i/100;
    DesiredLegLn = diag([1,multiplier,multiplier,1,1,multiplier])*InLegLn;
    
    error = (DesiredLegLn - ThLegLn);
    dlta = ThJac\error;
    ThPlatPose = NewOri(dlta, ThPlatPose);
    
    
    ThLegLn = zeros(6,1);
    ThJac = zeros(6,6);
    
    for ji = 1:6
        PlatPt = ThPlatPose*[legs(1:3,ji);1];          r = PlatPt(1:3) - ThPlatPose(1:3,4);
        BasePt = legs(4:6,ji);
        L = PlatPt(1) - BasePt(1);  M = PlatPt(2) - BasePt(2);  N = PlatPt(3) - BasePt(3) ;
        ThLegLn(ji) = sqrt(L^2+M^2+N^2);
        ThJac(ji,:) = [L,M,N,cross(r,[L;M;N])']/ThLegLn(ji);
    end
    
    if rank(ThJac)<6 ; disp('Platform configuration turned out to be singular'); return; end
    
end
