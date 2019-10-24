
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%Planar Event Reconstruction Algorithm%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script estimates the LRFs (Optical Model) for several Insert clinical modules by using a flood dataset (more that 100K events are 
% required) for each. The flood dataset is a .mat file: a N x M matlab matrix where:
% N = number of events
% M = number of acquisition channels 

% [ The dataset is the 'Frame' variable produced by 'INSERT_reorder' script. This script has to be run before 'LRFs_estimation' and
% the variable 'Frame' generated by the code has to be saved in 'Database'
% folder of 'PERA_PlanarReconstructionAlgorithm'.]

% LRFs that are produced by this script for each module are automatically
% saved in 'LRFs' folder.

close all
clear all
clc

%% CALIBRATION DATASETS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the name of the calibration dataset (saved in 'Database') to be used for  the calibration of each module:

Name{1} = 'UniformFlood_Nd01'; % Calibration dataset for module #1
Name{2} = 'UniformFlood_Nd02';
Name{3} = 'UniformFlood_Nd03';
Name{4} = 'UniformFlood_Nd04';
Name{5} = 'UniformFlood_Nd05';
Name{6} = 'UniformFlood_Nd06';
Name{7} = 'UniformFlood_Nd07';
Name{8} = 'UniformFlood_Nd08';
Name{9} = 'UniformFlood_Nd09';
Name{10} = 'UniformFlood_Nd10';
Name{11} = 'UniformFlood_Nd11';
Name{12} = 'UniformFlood_Nd12';
Name{13} = 'UniformFlood_Nd13';
Name{14} = 'UniformFlood_Nd14';
Name{15} = 'UniformFlood_Nd15';
Name{16} = 'UniformFlood_Nd16';
Name{17} = 'UniformFlood_Nd17';
Name{18} = 'UniformFlood_Nd18';
Name{19} = 'UniformFlood_Nd19';
Name{20} = 'UniformFlood_Nd20'; % Calibration dataset for module #20

% disp('%%%%%%%%%%%%%%%%%%%% PERA %%%%%%%%%%%%%%%%%%')
% disp('>>> CALIBRATION DATASETS <<<');
% n_modules = input('Number of modules to be calibrated: ');
% clc
% 
% for i = 1:n_modules
%     disp('>>> CALIBRATION DATASETS <<<');
%     disp(strcat('Write the name of calibration dataset file (.mat) for module ', num2str(i),' [e.g. Flood_Clinical ]:'));
%      Name{i} = input('', 's');
%      clc
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% image_show hast to be set to 1 to show the images produced by the Optical_Model_Estimator function
image_show = 1;


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

if overwrite_parameters
    load('Optical_model_parameters.mat') %gives 'Tune' struct variable as output
    load('Filt') 
else
    % if not defined by the user, the baseline is set as the mean value of signal from the detectors that collected less light');
    disp('>>> LRFs MANUAL SETTINGS <<<')
    Filt.baseline = input('Set the baseline value [ADC_channels] to be subtracted to data for modified centroid method: ');
    clc
    disp('>>> LRFs MANUAL SETTINGS <<<')
    Tune.dirty_noise = input('Set the arbitrary blurring noise on data [# of pixels] applied on reconstruction coordinates to facilitate LRFs convergence to optimal result: ');
    clc
    disp('>>> LRFs MANUAL SETTINGS <<<')
    Tune.iteration_exit = input('Set the number of iterations for the optical model estimation: ') ;
    clc
    disp('>>> LRFs MANUAL SETTINGS <<<')
    Tune.Num_rec = input('Set the number of event to be used for optical model estimation:');
    clc
end

load('LRFs_sampling')
Tune.pixel_LRF = LRFs_sampling.sampling;

addpath(strcat(pwd,'\Database'));
addpath(strcat(pwd,'\Functions'));
addpath(strcat(pwd,'\Models and Corrections'))

%% OPTICAL MODEL ESTIMATION (all modules)
for m = 2:size(Name,2)
    
	clearvars -except m Name Filt Tune LRFs_sampling image_show
    close all
    Dataset_filename = [Name{m}];
    

%% LOADING CALIBRATION DATASET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Frame = CalibrationDatasetLoad(Dataset_filename);


%% ENERGY SPECTRUM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the EnergySpectrum function computes the histogram of the sum of the
% signal from all channels for all the events of the dataset and plots it 
% (this is an uncalibrated energy spectrum)

output.energy = EnergySpectrum(Frame, 1); 
disp('>>> SELECTION of PHOTOPEAK EVENTS <<<')
disp('Select energy values belonging to the photopeak to extract events for LRFs estimation:')
disp('place the pointer on the lower limit of the range and click. Repeat the same for the upper limit.')

%%
% Auto Peak
pmod = mode(output.energy);
energyWindowed = [output.energy(output.energy > (pmod-0.1*pmod));...
    output.energy(output.energy < (pmod+0.1*pmod))];

[c,d] = hist(energyWindowed,floor(sqrt(size(Frame,1))));

fitting=GaussFit(d,c,pmod);
x_peak=round(fitting.b1); 
figure(1), hold on, plot(d,c,'xr'), hold off;

en_window_perc_width=0.1; 
% Select the energy range for data energy windowing for image filtering)
%-energy windowing
Filt.E_min=round((1-en_window_perc_width).*x_peak);
Filt.E_max=round((1+en_window_perc_width).*x_peak);

%% Manual Peak 
% [px,py]=ginput(2);
% close all
% clc
% % Select the energy range for data energy windowing for image filtering)
% Filt.E_min= min(px);%channels
% Filt.E_max= max(px);%channels


%% OPTICAL MODEL ESTIMATION: LRFs 
% LRFs estimation function
[LRFs, output] = Optical_Model_Estimator(Frame, LRFs_sampling.sampling, Tune, Filt, output, image_show);


%% SAVE LRFs and RECONSTRUCTION
% LRFs savepath
LRFs_savepath = strcat(pwd,'\LRFs\LRFs_',Dataset_filename,'.mat');
% Reconstruction savepath
output_savepath = strcat(pwd,'\Database_Reconstructions\rec_20191015_',Dataset_filename,'.mat');

save(LRFs_savepath,'LRFs');
save(output_savepath, 'output');

    clc
    disp(['end module ',num2str(m)])
    
end