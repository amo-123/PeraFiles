function [LRF] = LRF_Load_Spline(Num_channel,LRFs,Par)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
LRF=zeros(Num_channel,length(Par.y_knot),length(Par.x_knot));
n_div=7;
factor=1;
dim_x=size(Par.x_knot,2)/factor;
dim_y=size(Par.y_knot,2)/factor;
n_breaks_y = n_div;
n_breaks_x = 2*n_div;
ky = 3;
kx = 3;
gain=0;%pixel%%%%%%%%%%%REMOVE
knotsx = augknt(linspace(1,dim_x,n_breaks_x+1),kx);
knotsy = augknt(linspace(1,dim_y,n_breaks_y+1),ky);


x_bin = 1:1:dim_x; y_bin = 1:1:dim_y;

for i=1:Num_channel
    f=LRFs{i};
    values = spcol(knotsx,kx,x_bin)*f.coefs*spcol(knotsy,ky,y_bin).'; 
    values=values';
%     figure
% surf(values,'Linestyle','none');
% drawnow
    LRF(i,:,:)=reshape(values,1,length(Par.y_knot),length(Par.x_knot));
    
end
%to avoid negative or null values
% LRF(LRF<=1e-6)=1e-6;
end
