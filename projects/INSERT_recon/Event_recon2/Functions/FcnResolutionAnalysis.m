function [ FWHM , pos] = FcnResolutionAnalysis( centroid_figure,width,pixel,n_points )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

%% Dimensione foro collimatore
dim_coll=0.5;%mm

%% Calcolo della risoluzione intrinseca per una griglia
figure
B=centroid_figure;
imagesc(B);
colormap('hot')
axis equal

%Choose manually the point in the image
[x1,y1]=ginput(1);
close
%Point magnification

X_proj=B(ceil(y1-width*pixel):floor(y1+width*pixel),:);
Y_proj=B(:,ceil(x1-width*pixel):floor(x1+width*pixel));
% numero_eventi_X=sum(sum(X_proj));
% numero_eventi_Y=sum(sum(Y_proj));

figure
imagesc(X_proj);
colormap('hot')
axis equal
figure
imagesc(Y_proj);
colormap('hot')
axis equal
%Choose the very center of the point

line_x=sum(X_proj,1);
line_y=sum(Y_proj,2)';
x_coords=pixel*[1:length(line_x)]';
y_coords=pixel*[1:length(line_y)]';
linearity_correction = 0;
if linearity_correction == 1
    fitx=load(strcat(pwd,'\Linearity\interpolant_x.mat'))
    fity=load(strcat(pwd,'\Linearity\interpolant_y.mat'))
    x_coords = fitx.fitresult(x_coords);
y_coords = fity.fitresult(y_coords);
end


figure
plot(x_coords,smooth(line_x,5),'k','Linewidth',3)
axis([min(x_coords) max(x_coords) 0 1.1*max(smooth(line_y,3))])
xlabel('X [mm]','Fontname','Georgia','Fontsize',14,'Fontweight','bold')
ylabel('counts [a.u.]','Fontname','Georgia','Fontsize',14,'Fontweight','bold')
set(gca,'Fontname','Georgia','Fontsize',14,'Fontweight','bold')
x_lim=ginput(2);
x_lim=x_lim(:,1);

% close
figure
plot(y_coords,smooth(line_y,3),'k','Linewidth',3)
axis([min(y_coords) max(y_coords) 0 1.1*max(smooth(line_y,3))])
xlabel('Y [mm]','Fontname','Georgia','Fontsize',14,'Fontweight','bold')
ylabel('counts [a.u.]','Fontname','Georgia','Fontsize',14,'Fontweight','bold')
set(gca,'Fontname','Georgia','Fontsize',14,'Fontweight','bold')
y_lim=ginput(2);
y_lim=y_lim(:,1);

% close
range_x=(x_coords>=x_lim(1)).*(x_coords<=x_lim(2));
range_y=(y_coords>=y_lim(1)).*(y_coords<=y_lim(2));
fit_name=['gauss',num2str(n_points)];

fit_Xprojection=fit(x_coords(range_x==1),smooth(line_x(range_x==1)',3),fit_name);
fit_Yprojection=fit(y_coords(range_y==1),smooth(line_y(range_y==1)',3),fit_name);
centers_x=[];
centers_y=[];
std_x=[];
std_y=[];
for i=1:n_points
eval(['centers_x=[centers_x fit_Xprojection.b',num2str(i),'];'])
eval(['std_x=[std_x fit_Xprojection.c',num2str(i),'];'])
eval(['centers_y=[centers_y fit_Yprojection.b',num2str(i),'];'])
eval(['std_y=[std_y fit_Yprojection.c',num2str(i),'];'])
end

[centers_x,order_x]=sort(centers_x);
std_x=std_x(order_x);
[centers_y,order_y]=sort(centers_y);
std_y=std_y(order_y);
% std_x(1)=[];std_x(end)=[];
% std_y(1)=[];std_y(end)=[];

pos.x = centers_x;
pos.y = centers_y;
FWHM.x.ext=2.35482*std_x/sqrt(2);
FWHM.y.ext=2.35482*std_y/sqrt(2);
FWHM.x.int=sqrt(FWHM.x.ext.^2-dim_coll^2);
FWHM.y.int=sqrt(FWHM.y.ext.^2-dim_coll^2);

figure
plot(x_coords,smooth(line_x,3),'Linewidth',2)
hold on
plot(fit_Xprojection)

figure
plot(y_coords,smooth(line_y,3),'Linewidth',2)
hold on
plot(fit_Yprojection)

end

