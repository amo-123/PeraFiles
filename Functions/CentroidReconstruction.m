function [ x_rec,y_rec,Counts,energy_window ,number_det_ch] = CentroidReconstruction( Frame,Par,Filt,Energy,depict)

disp('Centroid Recontruction')
%Define number of pixels in X and Y for image sampling
pixel = 0.2;%mm
biny = ceil(Par.cryst_lung_y/pixel); %ceil allows to round the values of the elements of the vector between brackets to the nearest integers towards infinity
binx = ceil(Par.cryst_lung_x/pixel); 
%% Definizione variabili
%Copy energy filter values
E_min=Filt.E_min; %used to select only those events whose associated signal belongs to the gamma photopeak (in this way events at lower energies, associated to scattering/escape, are discarded)
E_max=Filt.E_max;
%Baseline definition
if isfield(Filt,'baseline')
    Baseline = Filt.baseline;
else
    Baseline = mean(mean(Frame));%channels
end
%Detection channels centers
X_mat=repmat(Par.X_c,size(Frame,1),1);%creates a matrix whose rows represent 
% the x-coordinates of centers of all the photodetectors according to the 
%selected convention => this is used in order to exploit matrix computations 
%in event position reconstruction
Y_mat=repmat(Par.Y_c,size(Frame,1),1);

%% Centroid method
%Baseline subraction
Frame_sub=Frame-Baseline;
Frame_sub(Frame_sub<0)=0; %Frame elements that take negative values after 
%baseline subtraction are set to 0

%Number of active detection channels per event, after baseline subtraction
number_det_ch=sum(Frame_sub>0,2);

%Centroid Formula : it uses the modified centroid method to reconstruct the
%coordinates for all the detected events. In order to select only those 
%events that belong to the gamma photopeak some filters are
%applied in order to discard "not acceptable" reconstructed events.
Energy_sub=sum(Frame_sub,2);
x_rec=sum((Frame_sub.*X_mat),2)./Energy_sub; %x_rec is the vector of the reconstructed x coordinates for ALL the events in the input dataset
y_rec=sum((Frame_sub.*Y_mat),2)./Energy_sub;

%% Filters
%energy window
F1=Energy>E_min;
F2=Energy<E_max;
%Coordinates reconstructed in the Field of View (FOV)
F_tempX = ((x_rec>=0).*(x_rec<=Par.cryst_lung_x))>0;
F_tempY = ((y_rec>=0).*(y_rec<=Par.cryst_lung_y))>0;
F_FOV = F_tempX.*F_tempY;
%Filter dependent on the number of detection channels activated (influences
%only the image)
F3=number_det_ch<=size(Frame,2);
F4=number_det_ch>=3;

%Image filter (for Centroid reconstruction image)
filter_image=find((F1.*F2.*F3.*F4.*F_FOV)==1);

y_image=y_rec(filter_image,:)';%x_image and y_image are variables used just for image display
x_image=x_rec(filter_image,:)';
%two "fake" events are added at the corners of the FOV to ensure the
%histogram to get all the FOV (they are not considered in the output)

y_image  = [y_image 0 Par.cryst_lung_y];
x_image  = [x_image 0 Par.cryst_lung_x];

[Counts,Centers]=hist3([y_image' x_image'],[biny binx]);
% Counts(1,1)=Counts(1,1)-1;
% Counts(end,end)=Counts(end,end)-1;

%% Image
if depict
    figure
    % immagine=flipud(immagine);
    imagesc([0 Par.cryst_lung_x],[0 Par.cryst_lung_y],Counts)
    set(gca,'Fontsize',15,'Fontname','TimesNewRoman','FontWeight','bold')
    xlabel('x [mm]')
    ylabel('y [mm]')
    hold on
    colormap hot
    axis equal
    axis tight
    colorb=colorbar;
    set(colorb,'Fontsize',12);
    set(colorb,'FontWeight','bold');
    set(colorb,'Fontname','Arial');
    title('Modified Centroid Reconstruction', 'Fontsize',14);
    drawnow;
end

%% MIA VERSIONE DEL PLOT 
% edges{1,2} = linspace(0, Par.cryst_lung_x,binx);
% edges{1,1} = linspace(0, Par.cryst_lung_y,biny);
% [Counts]= hist3([y_image' x_image'], 'Edges', edges);
% figure
% imagesc([0 Par.cryst_lung_x],[0 Par.cryst_lung_y],Counts);
% set(gca,'Fontsize',15,'Fontname','TimesNewRoman','FontWeight','bold')
% xlabel('x [mm]')
% ylabel('y [mm]')
% hold on
% colormap hot
% axis equal
% axis tight
% colorb=colorbar;
% set(colorb,'Fontsize',12);
% set(colorb,'FontWeight','bold');
% set(colorb,'Fontname','Arial');
% title('Modified Centroid Reconstruction', 'Fontsize',14);
% drawnow;
%     
%% Output reconstruction
%Filter energy window (for LRF estimation)
energy_window=find((F1.*F2.*F_FOV)==1);%It gives as output a vector containing 1 in those positions corresponsing to events that: a)match the energy filter b)whose positions are reconstructed within the camera FOV 
y_rec=y_rec(energy_window,:)';%x_rec and y_rec are reconstructed positions returned by this function 
x_rec=x_rec(energy_window,:)';

%% end
clc
end

