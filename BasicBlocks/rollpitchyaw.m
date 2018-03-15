clear variables;
close all;
clc;
%% This script find outs the Orientation matrix from roll pitch yaw convention
% sequence that is followed is, first yaw about x axis of base, followed by
% pitch about the y axis of base and finally a roll about z axis

syms alpha beta gamma real

thz = alpha;    % roll
thy = beta;     % pitch
thx = gamma;    % yaw

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