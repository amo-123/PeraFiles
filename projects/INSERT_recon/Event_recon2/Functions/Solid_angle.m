function [ solid_angle ] = Solid_angle( x,y,z, xi,yi,dx,dy )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% x=reshape(x,length(x),1);
% y=reshape(y,1,length(y));
SA1 = atan( (x-(xi-dx/2))*(y-(yi-dy/2)) / (z*sqrt((x -(xi-dx/2))^2+(y -(yi-dy/2))^2+z^2)) );
SA2 = atan( (x-(xi-dx/2))*(y-(yi+dy/2)) / (z*sqrt((x -(xi-dx/2))^2+(y -(yi+dy/2))^2+z^2)) );
SA3 = atan( (x-(xi+dx/2))*(y-(yi+dy/2)) / (z*sqrt((x -(xi+dx/2))^2+(y -(yi+dy/2))^2+z^2)) );
SA4 = atan( (x-(xi+dx/2))*(y-(yi-dy/2)) / (z*sqrt((x -(xi+dx/2))^2+(y -(yi-dy/2))^2+z^2)) );

solid_angle = SA1 - SA2 + SA3 - SA4; 
solid_angle = solid_angle/(4*pi);
end

