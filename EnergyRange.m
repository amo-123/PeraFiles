function [ enwind ] = EnergyRange(filename,filepath,L_msk)

%% Main to Split Nodes
%Code to load .data files (INSERT files) and split the acquisitions from different nodes into FRAME_NODES variable

% close all
% clear all
% clc

% if ~exist('current_path', 'var')
%     current_path = pwd;
% end
% current_path=strcat(current_path,'/');
addpath(strcat(pwd,'/MATLAB_modified_Main_Split'));
addpath(strcat(pwd,'/MATLAB_modified_Main_Split/Functions'));
addpath(strcat(pwd,'/MATLAB_modified_Main_Split/Geometries'));
addpath(strcat(pwd,'/MATLAB_modified_Main_Split/Functions/dataFunctions'));


%% Initializations

%determine if the user wants to open only the last num_events
num_events=5000000; %number of events to display (comment this line to see all the events)

%% Load

% FilterSpec = '*.data';
% [filename,filepath] = uigetfile([pwd,'/',FilterSpec], 'Select .data file', 'MultiSelect', 'off');

% [Datasize] = MS_getfilesize(filename,filepath,num_events);
% 
% AllData = cell(Datasize,20);
% 
% for ii = 1:Datasize
    %load function
    if exist('num_events','var')
        [Frame,Node,~,modality]=openDataFileMod(filename,filepath,num_events,1); % Change last argument for range of file
        %dispstrcat('last',{' '},num2str(num_events),{' '},'events loaded'))
    else
        [Frame,Node,~,modality]=openDataFileMod(filename,filepath);
        %disp'all events loaded')
    end
    
    %load proper reorder array and geometrical reorder
    [Phys_order]=INSERT_reorder(modality);
    Frame(:,Phys_order)=Frame;
    Node=round(Node);
    num_nodes=max(Node); %number of nodes
    
    %divide datasets of different nodes
%     FRAME_NODE{num_nodes,1}=0;
    FRAME_NODE = cell(num_nodes,1);
    
%     [~,ic1] = max( Frame, [], 2 ); 
%     Frame = Frame( ic1>6, : );

    for n=1:num_nodes
        try
        FRAME_NODE{n,1}=Frame(Node==n,:);
        catch 
        FRAME_NODE{n,1}= [];
        end
    end
    FRAME_NODE = FRAME_NODE(~cellfun('isempty',FRAME_NODE));
    num_nodes=length(FRAME_NODE);
    %dispstrcat('Modality:',32,modality,32,'-> Number of Nodes:',32,num2str(num_nodes)))

    %% Check number of events x node
    
    % for n=1:num_nodes
    %     Events_counts(n)=length(FRAME_NODE{n,1});
    % end
    % figure
    % plot(1:1:num_nodes,Events_counts,'-o','linewidth',2)
    % xlabel('#node')
    % ylabel('#events')
    % set(gca,'fontsize',15,'fontweight','bold')
    % grid on
    
    % Tot_events=sum(Events_counts);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MLEM  %%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
rmpath(strcat(pwd,'/MATLAB_modified_Main_Split'));
rmpath(strcat(pwd,'/MATLAB_modified_Main_Split/Geometries'));
rmpath(strcat(pwd,'/MATLAB_modified_Main_Split/Functions'));
rmpath(strcat(pwd,'/MATLAB_modified_Main_Split/Functions/dataFunctions'));


    % close all
    % clear all
    % clc
    
    %% DATASET
    % Set here the name of the DATASET to be reconstructed. The dataset should be organized in a N x M matrix, where N is the number of events
    % and M the number of detection channels. Such a dataset is created as 'Frame' variable by 'INSERT_reorder' script (his variable has to
    % be manually saved in 'Database' folder).
    
    enwind = zeros(20,2);
    
    for jj = 1:num_nodes
        % file_name = '20170403_module9_Co57_gain15_th30_HV35e4_flood_02';
        
        
        %% LRFs for STATISTICAL RECONSTRUCTION
        % Set here the name of the file containing the LRFs (optical model) produced by 'LRF_estimation' script for the considered module.
        % The LRFs file is automatically saved by 'LRFs_estimation' in "LRFs"
        % folder.
%         lrf_folder = './LRFs/';
%         lrf_files = dir(fullfile(lrf_folder,'*.mat'));
%         
%         LRFs_filename = lrf_files(jj);
        
        
        %% BASELINE for MODIFIED CENTROID RECONSTRUCTION
        % Set the baseline value is subtracted to the signals of all the channels.
        % This is an arbitrary value and has to be defined by observing the result obtained in terms of reconstructed image employing the
        % Centroid Method: the baseline value should be set so that reconstructed events cover the whole FOV.
        % Typical values: from 400 to 600 ADC_channels.
        
%        baseline = 600; % ADC_channels
        
        %% ENERGY CALIBRATION/RESOLUTION SETTING
        % Set 'En_resolution' to:
        % 1 --> show reconstructed energy spectrum, perform energy calibration and compute energy resolution
        % 0 --> just show reconstructed energy spectrum
%         En_resolution = 0;
        
        
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
        addpath(strcat(pwd,'/Geometries'))
        addpath(strcat(pwd,'/Functions'))
        addpath(strcat(pwd,'/Database'))
        load('INSERT_Clinical.mat','Par')% Geometrical Parameters are saved in "Par" struct
        
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
%         Filt.error_max = 1;
%         Filt.error_min = 0.0;
        
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
%         Filt.baseline= baseline;%channels
%         
%         % Choose the statistical method for reconstruction
%         Filt.Recon_Method='ML'; %Choose 'ML' (Maximum Likelihood) or 'LS' (Least Square Error) or 'WS' (Weighted Square Error)
        
        % Loading tuning parameters for the iterative method
        % call of all dataset in "Models and Corrections" (Tune parameters are
        % redundant in this script)
        addpath(strcat(pwd,'/Models_and_Corrections'))
        load('Optical_model_parameters.mat','Tune') %gives 'Tune' structured variable as output
        
        % load dataset
        % Frame = CalibrationDatasetLoad(file_name);
        Frame = FRAME_NODE{jj};
        
        [~,ic1] = max( Frame, [], 2 );
        Frame = Frame( ic1>6, : );
        %%
        %Neighbourhood mask (KE)
        %Baseline subraction
        Frame_sub=Frame-600;
        Frame_sub(Frame_sub<0)=0;
        
        if ( L_msk )
            nx = 6;  ny = 12;
            nv = size(Frame,1);
            [~,ii_max] = max(Frame,[],2);
            ix = mod( ii_max-1, nx ) + 1;
            iy = fix( (ii_max-1) / nx ) + 1;
            msk = zeros([nx,ny,nv]);
            ix1=ix-2; ix1(ix1<1)=1; ix2=ix+2; ix2(ix2>nx)=nx;
            iy1=iy-1; iy1(iy1<1)=1; iy2=iy+1; iy2(iy2>ny)=ny;
            jx1=ix-1; jx1(jx1<1)=1; jx2=ix+1; jx2(jx2>nx)=nx;
            jy1=iy-2; jy1(jy1<1)=1; jy2=iy+2; jy2(jy2>ny)=ny;
            for iv=1:nv
                msk( ix1(iv):ix2(iv), iy1(iv):iy2(iv), iv )=1;
                msk( jx1(iv):jx2(iv), jy1(iv):jy2(iv), iv )=1;
            end
            msk = reshape( msk, [nx*ny,nv] )';
            Frame_sub = Frame_sub .* msk;
        end
        
        %% ENERGY SPECTRUM
        % the EnergySpectrum function calculates the output energy (sum of the
        % signal from all channels) and plots the signal spectrum (uncalibrated energy spectrum
        output.energy1 = EnergySpectrum(Frame, 0);
        
        %disp'%%%%%%%%%%%%%%%%%%%% PERA %%%%%%%%%%%%%%%%%%')
        %disp'>>> SELECTION of EVENTS to be reconstructed by statistical method <<<')
        %disp'Select energy values related to events to be reconstructed by the statistical method:')
        %disp'place the pointer on the lower limit of the range and click. Repeat the same for the upper limit.')
        
        %%
        % Auto Peak
        pmod = mode(output.energy1);
        energyWindowed = [output.energy1(output.energy1 > (pmod-0.1*pmod));...
            output.energy1(output.energy1 < (pmod+0.1*pmod))];
        
        [c,d] = hist(energyWindowed,floor(sqrt(size(Frame,1))));
        [a,b] = hist(energyWindowed,floor(sqrt(size(Frame_sub,1))));
        
        fitting=GaussFit(d,c,pmod);
        x_peak=round(fitting.b1);
%         figure, plot(d,c,'xb');
%         figure, plot(b,a,'xr');
        en_window_perc_width=0.15;
        % Select the energy range for data energy windowing for image filtering)
        %-energy windowing
        E_min=round((1-en_window_perc_width).*x_peak);
        E_max=round((1+en_window_perc_width).*x_peak);
        
    enwind(jj,:) = [E_min,E_max];
        
    end
   
