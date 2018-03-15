function [ Pfnew ] = NewOri( dlta, Pfold )
% takes delta (dx,dx,dz,dthx,dthy,dthz) and Pfold (platform old position) 
% and returns updated position of platform Pfnew

th = norm(dlta(4:6));                       % amount of rotation

kBase = dlta(4:6)/th;                       % vector of direction cosines represented in base frame
kx = kBase(1); ky = kBase(2); kz = kBase(3);   % direction cosines
%%
% kNew = Pfold(1:3,1:3)*kBase;                
% kx = kNew(1); ky = kNew(2); kz = kNew(3);   % direction cosines
%%

c = cos(th); s = sin(th); v = 1-c;

R = [kx^2*v+c      kx*ky*v-kz*s    kx*kz*v+ky*s   ;
    kx*ky*v+kz*s   ky^2*v+c        ky*kz*v-kx*s   ;
    kx*kz*v-ky*s   ky*kz*v+kx*s    kz^2*v+c    ];

Rot = R*Pfold(1:3,1:3);

Pfnew = [Rot,           Pfold(1:3,4)+dlta(1:3);
        zeros(1,3),     1];
end