function [LRFs, output] = Optical_Model_Estimator(Frame, LRFs_sampling,Tune, Filt,output, show)
addpath(strcat(pwd,'\Geometries'))
addpath(strcat(pwd,'\Functions'))

%% LOAD CAMERA GEOMETRY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load gamma camera geometrical data from "Geometries" folder %%%%%%%%%%%%%
% Please change the .mat file name depending on the camera geometry %%%%%%%
% fundamental parameters to be set:
% Par.X_c := X-coordinates of the detectors centers
% Par.Y_c := Y-ccordinates of the detectors centers
% Par.cryst_lung_x := lenght of the camera scintillator in the X direction
% Par.cryst_lung_y := lenght of the camera scintillator in the Y direction
% num_detect_x := number of detector columns in the X direction
% num_detect_y := number of detector raws in the Y direction
load('INSERT_Clinical.mat')% Geometrical Parameters are saved in "Par" struct

%% SETTINGS for OPTICAL MODEL ESTIMATOR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define sampling parameters and arrays for the Optical Model Map
% The software manager suggest not to manouver this parameter. Higher
% values make the reconstruction faster but with a higher error of
% estimation. On the other side, low values of Par.sampling increase the
% reconstruction accuracy, but increase the overall computing time.
Par.sampling = LRFs_sampling;
Par=SamplingParameters(Par);

% Filter on the reconstruction error from the statistical method - only for
% visualization, NO effect on the optical model estimation! (this filter
% has no effect on the image reconstructed with the Centroid Method)
Filt.error_max = 1;
Filt.error_min = 0.0;

% Choose the statistical method for reconstruction
Filt.Recon_Method='ML'; %Choose 'ML' (Maximum Likelihood) or 'LS' (Least Square Error) or 'WS' (Weighted Square Error)

% LRFs sampling
Tune.pixel_LRF = Par.sampling;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% MODIFIED CENTROID RECONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the pixel dimension for the event image reconstructed with the
% Centroid Method
Par.pixel_image = 0.2;%mm (ARBITRARY PARAMETER)

% CentroidReconstruction plots the image reconstructed with the
% Centroid Method and returns the reconstructed coordinates of all the
% events that follow the filter conditions
[output.x_rec_CM,output.y_rec_CM,output.Centroid_Counts,Filt.energy_window] = CentroidReconstruction(Frame,Par,Filt,output.energy,show);


%% INITIALIZATION FOR OPTICAL MODEL ESTIMATION METHOD

%save initial dataset and work with a subset defined in the energy window
Frame_init = Frame;
Frame = Frame_init(Filt.energy_window,:);

% choose a random sub dataset from Frame for Optical model estimation made up of 
% Tune.Num_rec events and saves it into REC.Frame_rec. 
Filt.subset_index = randsample(size(Frame,1),Tune.Num_rec);
REC.Frame_rec=Frame(Filt.subset_index,:);

%Coordinates reconstructed from centroid method for the selected 
% set of events are given as input to the iterative method
REC.x_rec=output.x_rec_CM(Filt.subset_index);
REC.y_rec=output.y_rec_CM(Filt.subset_index);


%% ITERATIVE METHOD
end_condition = 1; 
internal_counter = 0;

while(end_condition)
    %% 1ST STEP: Optical Model (LRF) Fitting
    [LRFs] = OpticalModelEngine(REC,Par,Tune,1,1); 
    
    
    %% 2ND STEP: Sampling of the LRF function 
    [Par.LRF] = LRF_Load(size(Frame,2),LRFs,Par);


    %% 3RD STEP: Statistical Method for reconstruction
    [output.x_rec,output.y_rec,output.energy,output.error,Filt] = StatisticalMethod_for_LRF( Frame,Par,Tune,Filt );
    
    
    
    %% UPDATES: Update next iteration reconstruction coordinates and display the current reconstructed image
    REC.x_rec=output.x_rec;
    REC.y_rec=output.y_rec;
    REC.Frame_rec=Frame(Filt.subset_index,:);
    
    %% SHOW CENTROID RECONSTRUCTION
    if (show == 1)
         %Display Reconstruction Result
         Par.pixel = 0.2;%mm
         DisplayReconstruction(output,Par,Filt,Tune,1 );
    end
    
    
    %% CHECK EXIT CONDITION AND UPDATE INTERNAL COUNTER
    internal_counter = internal_counter +1;
    if internal_counter >= Tune.iteration_exit
        end_condition = 0;
    end
    
end

%% STATISTICAL RECONSTRUCTION and DISPLAY of all the events: 2D histogram of reconstructed positions for all the events
if (show == 1)
    Tune.Num_rec=size(Frame,1);
    [output.x_rec,output.y_rec,output.energy,output.error,Filt] = StatisticalMethod_for_LRF( Frame,Par,Tune,Filt );
    [output.Statistical_Counts]= DisplayReconstruction( output,Par,Filt,Tune,1 );
end

return
