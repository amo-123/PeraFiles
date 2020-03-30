function [ DY,DX] = fitting_max_ML( val_fit )
%FITTING_MAX_ML Summary of this function goes here
%   Detailed explanation goes here
Zmp=val_fit(1,1);
Z0p=val_fit(1,2);
Zpp=val_fit(1,3);
Z00=val_fit(2,2);
Zm0=val_fit(2,1);
Zp0=val_fit(2,3);
Zmm=val_fit(3,1);
Z0m=val_fit(3,2);
Zpm=val_fit(3,3);
a0=Z00;
a1=-0.5*Zm0+0.5*Zp0;
a2=-0.5*Z0m+0.5*Z0p;
a3=0.25*Zmm-0.25*Zmp-0.25*Zpm+0.25*Zpp;
a4=0.5*Zm0-Z00+0.5*Zp0;
a5=0.5*Z0m-Z00+0.5*Z0p;
a6=-0.25*Zmm+0.25*Zmp+0.5*Z0m-0.5*Z0p-0.25*Zpm+0.25*Zpp;
a7=-0.25*Zmm+0.5*Zm0-0.25*Zmp+0.25*Zpm-0.5*Zp0+0.25*Zpp;
a8=0.25*Zmm-0.5*Zm0+0.25*Zmp-0.5*Z0m+Z00-0.5*Z0p+0.25*Zpm-0.5*Zp0+0.25*Zpp;
%parte iterativa
DX=-a1/a4/2;
DY=-a2/a5/2;
for i=1:2% DX e DY si misurano in bin
 DX=-(a1+a3*DY+a7*DY^2)/(a4+a6*DY+a8*DY^2)/2;
 DY=-(a2+a3*DX+a6*DX^2)/(a5+a7*DX+a8*DX^2)/2;
end

end

