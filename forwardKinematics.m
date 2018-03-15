close all;
clear all;
% clc;

%% This script computes the forward kinmatics of stewart platform 
%   using newton raphson method for error control, screw based kinematic
%   jacobian and inverse kinematics

%%
% InPfOri = [Roll, Pitch, Yaw, X, Y, Z] initial platform
% position

InPfOri = [degtorad(0); degtorad(20); degtorad(20); 0; 0; 300];

% Each column of the matrix legs corresponds to a leg of stewart platform and 
% have the structure [px; py; pz; bx, by, bz]. p stands for platform and b 
% stands for base. (px, py, pz) are the coordinates of the leg attachment point
% on platform, represented in the frame rigidly connected to the center of top 
% platform. Similarly (bx, by, bz) are the coordinated of point on leg attached 
% to base, represented in the frame rigidly attached to the base centroid.

legs = [100,    50,     -87,    -100,   -50,    87;
        0,      87,     50,     0,      -87,    -50;
        0,      0,      0,      0,      0,      0;
        150,    75,     -130,   -150,   -75,    130;
        0,      130,    75,     0,      -130,   -75;
        0,      0,      0,      0,      0,      0];
%%
% Computing transformation matrix representing the top platform orientation wrt 
% base. Following the convention yaw pitch roll for orientation
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
% (P, Q, R), (P*, Q*, R*) and (L, M, N) notations are from the screw theory


% r is the vector joining platform origin to PlatPt represented in
% base frame
for i = 1:6
    % PlatPt is coordinate of leg end attached to platform, represented
    % in base frame
    PlatPt = InPlatPose*[legs(1:3,i);1];          
    BasePt = legs(4:6,i);

    % r is the vector from the base of the platform to leg attachment point on
    % the platform represented in base frame
    r = PlatPt(1:3) - InPlatPose(1:3,4);
    
    L = PlatPt(1) - BasePt(1);  M = PlatPt(2) - BasePt(2);  N = PlatPt(3) - BasePt(3) ;
    InLegLn(i) = sqrt(L^2+M^2+N^2);
    InJac(i,:) = [L,M,N,cross(r,[L;M;N])']/InLegLn(i);
end

if rank(InJac)<6 ; disp('Intial configuration of platform was singular'); return; end

clearvars -except legs InLegLn InPlatPose

%%

GuessPlatPose = InPlatPose;
%%
disp(GuessPlatPose)
disp(InLegLn)

for multiplier = 1.01:0.01:1.05

    % AcLegLn = Actual Leg Lengths [L1,L2,L3,L4,L5,L6]
    AcLegLn = diag([1,multiplier,multiplier,1,1,multiplier])*InLegLn;

    % newton raphson method 
    % guessing the pose, then computing the leglengths corresponding to the
    % guesses pose, then updating the guessed pose and so on till the
    % difference between the actual leg lengths and leg length corresponding to
    % the updated guessed pose is small. If the loop converges, the guessed
    % pose will give the forward kinematics corresponding to the actualleg
    % lengths 
    for i = 1:100                       
        ThLegLn = zeros(6,1);
        jac = zeros(6,6);               % with screws stacked horizontally   
        for ji = 1:6
            PlatPt = GuessPlatPose*[legs(1:3,ji);1];          r = PlatPt(1:3) - GuessPlatPose(1:3,4);
            BasePt = legs(4:6,ji);
            L = PlatPt(1) - BasePt(1);  M = PlatPt(2) - BasePt(2);  N = PlatPt(3) - BasePt(3) ;
            ThLegLn(ji) = sqrt(L^2+M^2+N^2);
            jac(ji,:) = [L,M,N,cross(r,[L;M;N])']/ThLegLn(ji);
        end

        error = AcLegLn - ThLegLn;
        if norm(error) < 1e-4
            break
            % breaking the loop in case the error in leg lengths is extremely
            % small
        end

        if rank(jac) < 6
            disp('Jacobian turned out to be singular....  Stopping')
            break
        end
        
        dlta = jac\error;
        GuessPlatPose = NewOri(dlta, GuessPlatPose);

    end
    disp(GuessPlatPose)
    disp(ThLegLn)
end
