%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%Planar Event Reconstruction Algorithm%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script performs event reconstruction of a dataset by using modified centroid algorithm and a statistical reconstruction algorithm.
% The script requires as inputs:
% 1) name of the dataset file
% 2) name of the LRFs (optical model) to be used for statistical reconstruction
% 3) data baseline value

% Reconstruction is automatically saved in 'Database_Reconstruction' folder. It contains:
% >> reconstructed coordinates by centroid algorithm
% >> reconstructed coordinates by statistical algorithm
% >> reconstructed energy by statistical algorithm
% >> reconstruction error by statistical algorithm
tic;
close all
clear all
clc

%% DATASET 
% Set here the name of the DATASET to be reconstructed. The dataset should be organized in a N x M matrix, where N is the number of events 
% and M the number of detection channels. Such a dataset is created as 'Frame' variable by 'INSERT_reorder' script (his variable has to 
% be manually saved in 'Database' folder). 

%file_name = '20170403_module9_Co57_gain15_th30_HV35e4_flood_02';
folder = uigetdir;
files = dir(fullfile(folder,'*.mat'));

LRFfolder = uigetdir;
LRFfiles = dir(fullfile(LRFfolder,'*.mat'));
for i = 1:length(files) %change
fn = ['\London\X\',files(i).name];
file_name = fn(1:end-4);
%filepath = [files(i).folder,'\'];

%file_name = 'UniformFlood_Nd01';

%% LRFs for STATISTICAL RECONSTRUCTION
% Set here the name of the file containing the LRFs (optical model) produced by 'LRF_estimation' script for the considered module.
% The LRFs file is automatically saved by 'LRFs_estimation' in "LRFs"
% folder.

%LRFs_filename = 'LRF_2017-10-26_bulmaraw_H01_F';

LRFs_filename = LRFfiles(i).name(1:end-4);
%LRFs_filename = 'LRFs_Node20_TestFlood';
%% BASELINE for MODIFIED CENTROID RECONSTRUCTION
% Set the baseline value is subtracted to the signals of all the channels. 
% This is an arbitrary value and has to be defined by observing the result obtained in terms of reconstructed image employing the 
% Centroid Method: the baseline value should be set so that reconstructed events cover the whole FOV. 
% Typical values: from 400 to 600 ADC_channels.

baseline = 500; % ADC_channels

%% ENERGY CALIBRATION/RESOLUTION SETTING
% Set 'En_resolution' to:
% 1 --> show reconstructed energy spectrum, perform energy calibration and compute energy resolution
% 0 --> just show reconstructed energy spectrum
En_resolution = 0;


%% PARAMETER SETTING for RECONSTRUCTION
% load gamma camera geometrical data from "Geometries" folder
% Please change the .mat file name depending on the camera geometry %%%%%%%
% fundamental parameters to be set:
% Par.X_c := X-coordinates of the detectors centers
% Par.Y_c := Y-ccordinates of the detectors centers
% Par.cryst_lung_x := lenght of the camera scintillator in the X direction
% Par.cryst_lung_y := lenght of the camera scintillator in the Y direction
% num_detect_x := number of detector columns in the X direction
% num_detect_y := number of detector raws in the Y direction
addpath(strcat(pwd,'\Geometries'))
addpath(strcat(pwd,'\Functions'))
addpath(strcat(pwd,'\Database'))
load('INSERT_Clinical.mat')% Geometrical Parameters are saved in "Par" struct

% Define sampling parameters and arrays for the Optical Model Map
% The software manager suggest not to manouver this parameter. Higher
% values make the reconstruction faster but with a higher error of
% estimation. On the other side, low values of Par.sampling increase the
% reconstruction accuracy, but increase the overall computing time.
Par.sampling=0.2;%mm 
Par=SamplingParameters(Par);

% Pixel dimension for the event image reconstructed with the
% Centroid Method
Par.pixel_image = 0.2;%mm (ARBITRARY PARAMETER)

% Filter on the reconstruction error from the statistical method - only for
% visualization, NO effect on the optical model estimation! (this filter
% has no effect on the image reconstructed with the Centroid Method)
Filt.error_max = 1; 
Filt.error_min = 0.0;

%%SPATIAL FILTERS : uncomment the following code to show reconstructed events belonging to the
%%UFOV
% UFOV_x = 90; %mm
% UFOV_y = 40; %mm 
% Filt.x_rec_min = (Par.cryst_lung_x - UFOV_x)/2; %mm
% Filt.x_rec_max = Par.cryst_lung_x - (Par.cryst_lung_x - UFOV_x)/2; %mm
% Filt.y_rec_min = (Par.cryst_lung_y - UFOV_y)/2; %mm
% Filt.y_rec_max = Par.cryst_lung_y - (Par.cryst_lung_y - UFOV_y)/2; %mm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER INITIALIZATION 
% Baseline subtraction for modified centroid method
%NB: if the following line is commented, the algorithm uses as baseline
%value the mean value of measured signals
Filt.baseline= baseline;%channels

% Choose the statistical method for reconstruction
Filt.Recon_Method='ML'; %Choose 'ML' (Maximum Likelihood) or 'LS' (Least Square Error) or 'WS' (Weighted Square Error)

% Loading tuning parameters for the iterative method
% call of all dataset in "Models and Corrections" (Tune parameters are
% redundant in this script)
addpath(strcat(pwd,'\Models and Corrections'))
load('Optical_model_parameters.mat') %gives 'Tune' structured variable as output

% load dataset
Frame = CalibrationDatasetLoad(file_name);
    [~,ic1] = max( Frame, [], 2 ); 
    Frame = Frame( ic1>6, : );

% Frame(:,14) = 0;

%% ENERGY SPECTRUM 
% the EnergySpectrum function calculates the output energy (sum of the
% signal from all channels) and plots the signal spectrum (uncalibrated energy spectrum
output.energy1 = EnergySpectrum(Frame, 1); 

disp('%%%%%%%%%%%%%%%%%%%% PERA %%%%%%%%%%%%%%%%%%')
disp('>>> SELECTION of EVENTS to be reconstructed by statistical method <<<')
disp('Select energy values related to events to be reconstructed by the statistical method:')
disp('place the pointer on the lower limit of the range and click. Repeat the same for the upper limit.')
%%
% Auto Peak
pmod = mode(output.energy1);
energyWindowed = [output.energy1(output.energy1 > (pmod-0.1*pmod));...
    output.energy1(output.energy1 < (pmod+0.1*pmod))];

[c,d] = hist(energyWindowed,floor(sqrt(size(Frame,1))));

fitting=GaussFit(d,c,pmod);
x_peak=round(fitting.b1); 
figure, plot(d,c,'xb');

en_window_perc_width=0.05; 
% Select the energy range for data energy windowing for image filtering)
%-energy windowing
Filt.E_min=round((1-en_window_perc_width).*x_peak);
Filt.E_max=round((1+en_window_perc_width).*x_peak);
%% Manual peak
% [px,py]=ginput(2);
% close all
% clc
% % % Select the energy range for data energy windowing for image filtering)
%  Filt.E_min= min(px);%channels
%  Filt.E_max= max(px);%channels


%% CENTROID RECONSTRUCTION 

% CentroidReconstruction plots the image reconstructed with the
% Centroid Method and returns the reconstructed coordinates of all the
% events that follow the filter conditions
[output.x_rec_CM,output.y_rec_CM,output.Centroid_Counts,Filt.energy_window] = CentroidReconstruction(Frame,Par,Filt,output.energy1,1);

% Save initial dataset and work with a subset defined by the selected
% energy window. Comment the two lines below to work with the wholw
% dataset.
Frame_init = Frame;
Frame = Frame_init(Filt.energy_window,:);


%% LOAD LRFs (Optical Model)

% Loading and sampling of the analytical selected LRF
load(strcat(pwd,'\LRFs\', LRFs_filename,'.mat'));
[Par.LRF] = LRF_Load(size(Frame,2),LRFs,Par);


%% STATISTICAL METHOD RECONSTRUCTION 

%overwrite Tune.Num_rec variable in order to process all the events instead
%of a subset
Tune.Num_rec=size(Frame,1);

% Employment of the statistical method: it is used to reconstruct the
% coordinates of the filtered events (according to the Energy filter)
[output.x_rec,output.y_rec,output.energy,output.error,Filt] = StatisticalMethod(Frame,Par,Tune,Filt );


%% HISTOGRAM of the RECONSTRUCTION ERROR of STATISTICAL RECONSTRUCTION
figure
hold on
set(gca,'Fontsize',14,'Fontname','Arial','FontWeight','bold')
title('Reconstruction error (Statistical Method)','Fontname','Arial','FontWeight','bold')
hist(output.error,100);


%% DISPLAY STATISTICAL RECONSTRUCTION
% 2D histogram of the reconstructed positions (XY)

Par.pixel = 0.2;%mm
[output.Statistical_Counts]=DisplayReconstruction(output,Par,Filt,Tune,0,0);
hold on
set(gca,'Fontsize',14,'Fontname','Arial','FontWeight','bold')
title('Statistical Reconstruction','Fontname','Arial','FontWeight','bold')


%% RECONSTRUCTED ENERGY HISTOGRAM, ENERGY CALIBRATION and RESOLUTION
energy_max = 1.5e+5;%ADC channels; valore al quale si vuole venga tagliato il plot dello spettro di energia

if En_resolution == 1
    %Energy peaks for spectrum calibration (these values have to be set in
    %order to compute energy resolution)
    disp('>>> ENERGY CALIBRATION AND RESOLUTION <<<')
    Filt.E_calibration_min = input('Energy of lower energy peak [keV]: ');%keV
    Filt.E_calibration_max =  input('Energy of higher energy peak [keV]: '); %keV
    clc
    
end

[Energy_Resolution, fitresult_energy_calibration, fit_ADC] = Reconstructed_Energy(output,energy_max,Filt,0,1,En_resolution);

if En_resolution == 1
       calibration_line = coeffvalues(fitresult_energy_calibration);

    disp(strcat('Energy calibration line is: Energy[keV] = ', num2str(calibration_line(1)),' * Energy[ADC_channels] + (', num2str(calibration_line(2)), ')'))
    disp(strcat('Energy resolution (FWHM of reconstructed energy spectrum) is: ', num2str(Energy_Resolution), '%'))

end

%% SAVE RECONSTRUCTION
% Reconstructed event positions (X,Y), energy, reconstruction error are
% stored in the 'output' structured variable. This one is saved in
% Database_Reconstructions folder

% Reconstruction savepath
output_savepath = strcat(pwd,'\Database_Reconstructions\Rec_',file_name(11:end),'.mat');
% Save reconstruction
save(output_savepath, 'output');
end
toc;