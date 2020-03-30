
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%Planar Event Reconstruction Algorithm%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clean(); %clears all variables and closes all open figures
clear


%% CALIBRATION DATASET NAME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%flood dataset of more than 100 Kevents, organized in a N x M matrix, where N is the number of events and M the
% number of detection channels
Dataset_filename = '20170403_module9_Co57_gain15_th30_HV35e4_flood_02';

%% OPTICAL MODEL ESTIMATOR SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overwrite_parameters = NaN;
while (overwrite_parameters ~= 0 && overwrite_parameters ~= 1)
    disp('%%%%%%%%%%%%%%%%%%%% PERA %%%%%%%%%%%%%%%%%%')
    disp('>>> LRFs SETTINGS <<<')
    disp('Type')
    disp('1: to use the default settings for optical model estimation')
    disp('0: to change the settings for optical model estimation manually')
    overwrite_parameters = input('');
    clc
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% image_show hast to be set to 1 to show the images produced by the Optical_Model_Estimator function
image_show = 0;

if overwrite_parameters
    load('Optical_model_parameters.mat') %gives 'Tune' struct variable as output
    load('Filt') 
else
    %if not defined by the user, the baseline is set as the mean value of signal from the detectors that collected less light');
    disp('2) LRFs MANUAL SETTINGS')
    Filt.baseline = input('Set the baseline value [ADC_channels] to be subtracted to data for modified centroid method: ');
    clc
    disp('2) LRFs MANUAL SETTINGS')
    Tune.dirty_noise = input('Set the arbitrary blurring noise on data [# of pixels] applied on reconstruction coordinates to facilitate LRFs convergence to optimal result: ');
    clc
    disp('2) LRFs MANUAL SETTINGS')
    Tune.iteration_exit = input('Set the number of iterations for the optical model estimation: ') ;
    clc
    disp('2) LRFs MANUAL SETTINGS')
    Tune.Num_rec = input('Set the number of event to be used for optical model estimation:');
    clc
end

load('LRFs_sampling')
Tune.pixel_LRF = LRFs_sampling.sampling;

%% LOADING CALIBRATION DATASET
addpath(strcat(pwd,'\Database'));
addpath(strcat(pwd,'\Functions'));
addpath(strcat(pwd,'\Models and Corrections'))
% load CALIBRATION DATASET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Frame = CalibrationDatasetLoad(Dataset_filename);

%% ENERGY SPECTRUM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the EnergySpectrum function computes the histogram of the sum of the
% signal from all channels for all the events of the dataset and plots it 
% (this is an uncalibrated energy spectrum)
disp(' >>> SELECTION of PHOTOPEAK EVENTS <<<')
output.energy = EnergySpectrum(Frame, 1); 
disp('Select energy values belonging to the photopeak to extract events for LRFs estimation:')
disp('place the pointer on the lower limit of the range and click. Repeat the same for the upper limit.')

[px,py]=ginput(2);
close all
clc
% Select the energy range for data energy windowing for image filtering)
Filt.E_min= min(px);%channels
Filt.E_max= max(px);%channels


%% OPTICAL MODEL ESTIMATION: LRFs 


% LRFs savepath
LRFs_savepath = strcat(pwd,'\LRFs\',Dataset_filename,'.mat');

% LRFs estimation function
[LRFs, output] = Optical_Model_Estimator(Frame, LRFs_sampling.sampling, Tune, Filt, output, image_show);

save(LRFs_savepath,'LRFs');
