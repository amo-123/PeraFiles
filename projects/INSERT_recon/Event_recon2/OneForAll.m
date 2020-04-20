function [ NodeData, AllData ] = OneForAll(filepath,filename,uflag,Ufilepath,Ufilename,EW,normflag,killchnl, outname)
% filepath: data file folder
% filename: data file name
% uflag: use flood spectra 
% Ufilepath: flood file folder
% Ufilename: flood file name 
% EW: Energy window file, if manual is required 
% killchnl: matrix of pk targets 
% outname: unique file identifier for output 

% FilterSpec = '*.data';

%[Ufilename,Ufilepath] = uigetfile([pwd,'\',FilterSpec], 'Select .data file', 'MultiSelect', 'off');
if uflag
 %   Ufilepath = '/SAN/inm/FDG/amoINSERT/Week_1/20190306/Flood/';
   
 %   Ufilename = '5ml1Mbq_Tc99m_flood_Tm10_hv35_gain12_th30_all_2min_00.data';
    [ enwind ] = EnergyRange(Ufilename,Ufilepath,1,killchnl);
    disp('EW found from U data');
elseif exist('EW','var') == 1
        load(EW,'enwind'); 
        disp('Manual EW');
else
    %load('EnergyWindows.mat','enwind');
   %load('/SAN/inm/FDG/amoINSERT/INSERT/PeraFiles/PeraFiles/EnergyWindows/Cylinder_02_EW.mat','enwind');
    %load('/SAN/inm/FDG/amoINSERT/INSERT/PeraFiles/PeraFiles/EnergyWindows/EW_Hoffman2D.mat','enwind');
%   load(EW,'enwind'); 
    [ enwind ] = EnergyRange(filename,filepath,1,killchnl);
    disp('EW found from data');
end

%% Main to Split Nodes
%Code to load .data files (INSERT files) and split the acquisitions from different nodes into FRAME_NODES variable

% close all
% clear all
% clc

% if ~exist('current_path', 'var')
%     current_path = pwd;
% end
% current_path=strcat(current_path,'\');
% <<<<<<< Updated upstream
% addpath('E:\TestLRF\PERA_PlanarReconstructionAlgorithm\PeraScripts\MATLAB_modified_Main_Split');
% addpath('E:\TestLRF\PERA_PlanarReconstructionAlgorithm\PeraScripts\MATLAB_modified_Main_Split\Functions');
% addpath('E:\TestLRF\PERA_PlanarReconstructionAlgorithm\PeraScripts\MATLAB_modified_Main_Split\Geometries');
% addpath('E:\TestLRF\PERA_PlanarReconstructionAlgorithm\PeraScripts\MATLAB_modified_Main_Split\Functions\dataFunctions');
% =======
addpath(strcat(pwd,'/MATLAB_modified_Main_Split'));
addpath(strcat(pwd,'/MATLAB_modified_Main_Split/Functions'));
addpath(strcat(pwd,'/MATLAB_modified_Main_Split/Geometries'));
addpath(strcat(pwd,'/MATLAB_modified_Main_Split/Functions/dataFunctions'));


%% Initializations

%determine if the user wants to open only the last num_events
num_events=5000000; %number of events to display (comment this line to see all the events)
n_chan = 72;
%% Load

%[filename,filepath] = uigetfile([pwd,'\',FilterSpec], 'Select .data file', 'MultiSelect', 'off');
% filepath = 'E:\Week 2\20190313_part2\CylinderPhantom\';
% filename = 'cylinder_50Mbq_Tc99m_Tm10_gain12_th30_hv35_2kill_06.data';

[Datasize] = MS_getfilesize(filename,filepath,num_events);

AllData = cell(Datasize-1,20);
AllFrame = cell(1,20);
X_rec = cell(1,20);
Y_rec = cell(1,20);

for ii = 1:Datasize - 1
    %load function
    if exist('num_events','var')
        [Frame,Node,~,modality]=openDataFileMod(filename,filepath,num_events,ii); % Change last argument for range of file
        disp(strcat('batch #', num2str(ii),' event loaded'));
    else
        [Frame,Node,~,modality]=openDataFile(filename,filepath);
        %disp('all events loaded')
    end
    
    %load proper reorder array and geometrical reorder
    [Phys_order]=INSERT_reorder(modality);
    Frame(:,Phys_order)=Frame;
    Node=round(Node);
    num_nodes=max(Node); %number of nodes
    
    %divide datasets of different nodes
%     FRAME_NODE{num_nodes,1}=0;
    FRAME_NODE = cell(num_nodes,1);

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
% <<<<<<< Updated upstream
%     rmpath('E:\TestLRF\PERA_PlanarReconstructionAlgorithm\PeraScripts\MATLAB_modified_Main_Split');
% rmpath('E:\TestLRF\PERA_PlanarReconstructionAlgorithm\PeraScripts\MATLAB_modified_Main_Split\Functions');
% rmpath('E:\TestLRF\PERA_PlanarReconstructionAlgorithm\PeraScripts\MATLAB_modified_Main_Split\Geometries');
% rmpath('E:\TestLRF\PERA_PlanarReconstructionAlgorithm\PeraScripts\MATLAB_modified_Main_Split\Functions\dataFunctions');
% =======
    rmpath(strcat(pwd,'/MATLAB_modified_Main_Split'));
rmpath(strcat(pwd,'/MATLAB_modified_Main_Split/Functions'));
rmpath(strcat(pwd,'/MATLAB_modified_Main_Split/Geometries'));
rmpath(strcat(pwd,'/MATLAB_modified_Main_Split/Functions/dataFunctions'));


    % close all
    % clear all
    % clc
    
    %% DATASET
    % Set here the name of the DATASET to be reconstructed. The dataset should be organized in a N x M matrix, where N is the number of events
    % and M the number of detection channels. Such a dataset is created as 'Frame' variable by 'INSERT_reorder' script (his variable has to
    % be manually saved in 'Database' folder).

    for jj = 1:num_nodes
        disp(strcat('Node #', num2str(jj),' Calculating...'));
        % file_name = '20170403_module9_Co57_gain15_th30_HV35e4_flood_02';

%         if length(FRAME_NODE) >= 13 %19
%             %FRAME_NODE{19,1}(:,14) = 0;
%             n_chan = 72;  ny = 6; % nz = 12;
%             jj_chan = (0:n_chan-1);
%             ii_chan = ( (mod(jj_chan,ny)>2) & (fix(jj_chan/ny)>5) );
%             FRAME_NODE{13,1}(:,ii_chan) = 0;
%         end
        for kk = 1:n_chan
            if killchnl(jj,kk) == 1 || killchnl(jj,kk) == 2
                FRAME_NODE{jj,1}(:,kk) = zeros(size(FRAME_NODE{jj,1}(:,kk)));
                disp(strcat('Node #', num2str(jj),' Channel #', num2str(kk), ' Killed.'));
            end
        end

        
        %% LRFs for STATISTICAL RECONSTRUCTION
        % Set here the name of the file containing the LRFs (optical model) produced by 'LRF_estimation' script for the considered module.
        % The LRFs file is automatically saved by 'LRFs_estimation' in "LRFs"
        % folder.
        lrf_folder = './LRFs/';
        lrf_files = dir(fullfile(lrf_folder,'*.mat'));
        if num_nodes == 1
            LRFs_filename = lrf_files(n).name(1:end-4);
        else 
            LRFs_filename = lrf_files(jj).name(1:end-4);
        end
        
        %% BASELINE for MODIFIED CENTROID RECONSTRUCTION
        % Set the baseline value is subtracted to the signals of all the channels.
        % This is an arbitrary value and has to be defined by observing the result obtained in terms of reconstructed image employing the
        % Centroid Method: the baseline value should be set so that reconstructed events cover the whole FOV.
        % Typical values: from 400 to 600 ADC_channels.
        
        baseline = 600; % ADC_channels
        
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
        load('INSERT_Clinical.mat','Par');% Geometrical Parameters are saved in "Par" struct
        
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
        Filt.L_msk = 1;
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

        addpath(strcat(pwd,'/Models_and_Corrections'));

        load('Optical_model_parameters.mat','Tune'); %gives 'Tune' structured variable as output
        
        % load dataset
        % Frame = CalibrationDatasetLoad(file_name);
        Frame = FRAME_NODE{jj};
        
        %% ENERGY SPECTRUM
        % the EnergySpectrum function calculates the output energy (sum of the
        % signal from all channels) and plots the signal spectrum (uncalibrated energy spectrum
        %output.energy1 = EnergySpectrum(Frame, 0);
        
        %disp('%%%%%%%%%%%%%%%%%%%% PERA %%%%%%%%%%%%%%%%%%')
        %disp('>>> SELECTION of EVENTS to be reconstructed by statistical method <<<')
        %disp('Select energy values related to events to be reconstructed by the statistical method:')
        %disp('place the pointer on the lower limit of the range and click. Repeat the same for the upper limit.')
        
        %%
%         % Auto Peak
%         pmod = mode(output.energy1);
%         energyWindowed = [output.energy1(output.energy1 > (pmod-0.1*pmod));...
%             output.energy1(output.energy1 < (pmod+0.1*pmod))];
%         
%         [c,d] = hist(energyWindowed,floor(sqrt(size(Frame,1))));
%         
%         fitting=GaussFit(d,c,pmod);
%         x_peak=round(fitting.b1);
%         %figure, plot(d,c,'xb');
%         
%         en_window_perc_width=0.05;
%         % Select the energy range for data energy windowing for image filtering)
%         %-energy windowing
%         Filt.E_min=round((1-en_window_perc_width).*x_peak);
%         Filt.E_max=round((1+en_window_perc_width).*x_peak);
        Filt.E_min = enwind(jj,1);
        Filt.E_max = enwind(jj,2);
        
        %% Maual Peak
        % [px,py]=ginput(2);
        % close all
        % clc
        % % Select the energy range for data energy windowing for image filtering)
        % Filt.E_min= min(px);%channels
        % Filt.E_max= max(px);%channels
        
        %% Normalise Channel Events
        disp('Channel Normalisation');
        if normflag == 1
            ch = sum(Frame);
            ch1 = ones(1,n_chan);
            if ( max(ch) > 0 )
                ii_chan = ( ch > 0 );
                ch1(ii_chan) = ch(ii_chan) / mean( ch(ii_chan) );
            end
            
            
            nrm = ones(size(Frame,1),1) * ch1;
            Frame = Frame ./ nrm;
        end
        
       
        
        %% CENTROID RECONSTRUCTION
        
        % CentroidReconstruction plots the image reconstructed with the
        % Centroid Method and returns the reconstructed coordinates of all the
        % events that follow the filter conditions
       % [output.x_rec_CM,output.y_rec_CM,output.Centroid_Counts,Filt.energy_window] = CentroidReconstruction(Frame,Par,Filt,output.energy1,0);
        %[output.x_rec_CM,output.y_rec_CM,output.Centroid_Counts,Filt.energy_window] = CentroidReconstructionKE(Frame,Par,Filt);
        disp('Centroid Recon')
        [output.x_rec_CM,output.y_rec_CM,output.Centroid_Counts,Filt.energy_window] = CentroidRecon(Frame,Par,Filt);

        
        % Save initial dataset and work with a subset defined by the selected
        % energy window. Comment the two lines below to work with the wholw
        % dataset.
        Frame_init = Frame;
        Frame = Frame_init(Filt.energy_window,:);
        
        
        %% LOAD LRFs (Optical Model)
        
        % Loading and sampling of the analytical selected LRF
        load(strcat(pwd,'/LRFs/', LRFs_filename,'.mat'),'LRFs');
        disp(strcat('Loading LRF -',LRFs_filename));
        
        for kk = 1:n_chan
            if killchnl(jj,kk) == 2
                LRFs{kk}.a = 0.02;
                LRFs{kk}.b = 16;
                LRFs{kk}.c = 16;
                disp(strcat('LRF Node #', num2str(jj),' Channel #', num2str(kk), ' Killed.'));
            end
        end
        
        [Par.LRF] = LRF_Load(size(Frame,2),LRFs,Par);
        

        
%         if normflag == 1
%             ch = sum(Frame);
%             ch1 = ones(1,n_chan);
%             if ( max(ch) > 0 )
%                 ii_chan = ( ch > 0 );
%                 ch1(ii_chan) = ch(ii_chan) / mean( ch(ii_chan) );
%             end
%             
%             
%             nrm = ones(size(Frame,1),1) * ch1;
%             Frame = Frame ./ nrm;
%         end
%%  DOI calculation 

        
        %% STATISTICAL METHOD RECONSTRUCTION
        
        %overwrite Tune.Num_rec variable in order to process all the events instead
        %of a subset
        Tune.Num_rec=size(Frame,1);
        
        % Employment of the statistical method: it is used to reconstruct the
        % coordinates of the filtered events (according to the Energy filter)
        [output.x_rec,output.y_rec,output.energy,output.error,Filt] = StatisticalMethod(Frame,Par,Tune,Filt );
        
        
        %% HISTOGRAM of the RECONSTRUCTION ERROR of STATISTICAL RECONSTRUCTION
%         figure
%         hold on
%         set(gca,'Fontsize',14,'Fontname','Arial','FontWeight','bold')
%         title('Reconstruction error (Statistical Method)','Fontname','Arial','FontWeight','bold')
%         hist(output.error,100);
        
        
        %% DISPLAY STATISTICAL RECONSTRUCTION
        % 2D histogram of the reconstructed positions (XY)
        
        Par.pixel = 0.2;%mm
        disp('ML Recon');
        [output.Statistical_Counts]=DisplayReconstruction(output,Par,Filt,Tune,0,0);
%         hold on
%         set(gca,'Fontsize',14,'Fontname','Arial','FontWeight','bold')
%         title('Statistical Reconstruction','Fontname','Arial','FontWeight','bold')
        
        
        %% RECONSTRUCTED ENERGY HISTOGRAM, ENERGY CALIBRATION and RESOLUTION
%         energy_max = 1.5e+5;%ADC channels; valore al quale si vuole venga tagliato il plot dello spettro di energia
%         
%         if En_resolution == 1
%             %Energy peaks for spectrum calibration (these values have to be set in
%             %order to compute energy resolution)
%             %disp'>>> ENERGY CALIBRATION AND RESOLUTION <<<')
%             Filt.E_calibration_min = input('Energy of lower energy peak [keV]: ');%keV
%             Filt.E_calibration_max =  input('Energy of higher energy peak [keV]: '); %keV
%             clc
%             
%         end
        
%         [Energy_Resolution, fitresult_energy_calibration, ~] = Reconstructed_Energy(output,energy_max,Filt,0,1,En_resolution);
%         
%         if En_resolution == 1
%             calibration_line = coeffvalues(fitresult_energy_calibration);
%             
%             %dispstrcat('Energy calibration line is: Energy[keV] = ', num2str(calibration_line(1)),' * Energy[ADC_channels] + (', num2str(calibration_line(2)), ')'))
%             %dispstrcat('Energy resolution (FWHM of reconstructed energy spectrum) is: ', num2str(Energy_Resolution), '%'))
%             
%         end
        
        %% SAVE RECONSTRUCTION
        % Reconstructed event positions (X,Y), energy, reconstruction error are
        % stored in the 'output' structured variable. This one is saved in
        % Database_Reconstructions folder
        AllData{ii,jj} = output;
        if ii == 1
            AllFrame{ii,jj} = Frame;
            X_rec{ii,jj} =  output.x_rec;
            Y_rec{ii,jj} =  output.y_rec;
        end
    end
% <<<<<<< Updated upstream
%             rmpath(strcat(pwd,'\Geometries'));
%         rmpath(strcat(pwd,'\Functions'));
%         rmpath(strcat(pwd,'\Database'));
%         
%     addpath('E:\TestLRF\PERA_PlanarReconstructionAlgorithm\PeraScripts\MATLAB_modified_Main_Split');
% addpath('E:\TestLRF\PERA_PlanarReconstructionAlgorithm\PeraScripts\MATLAB_modified_Main_Split\Functions');
% addpath('E:\TestLRF\PERA_PlanarReconstructionAlgorithm\PeraScripts\MATLAB_modified_Main_Split\Geometries');
% addpath('E:\TestLRF\PERA_PlanarReconstructionAlgorithm\PeraScripts\MATLAB_modified_Main_Split\Functions\dataFunctions');
% =======
            rmpath(strcat(pwd,'/Geometries'));
        rmpath(strcat(pwd,'/Functions'));
        rmpath(strcat(pwd,'/Database'));
        
addpath(strcat(pwd,'/MATLAB_modified_Main_Split'));
addpath(strcat(pwd,'/MATLAB_modified_Main_Split/Functions'));
addpath(strcat(pwd,'/MATLAB_modified_Main_Split/Geometries'));
addpath(strcat(pwd,'/MATLAB_modified_Main_Split/Functions/dataFunctions'));


end

% output_savepath = strcat(pwd,'\Database_Reconstructions\All_nodes_',filename(1:end-5),'.mat');
% save(output_savepath,'AllData')


% Reconstruction savepath
% output_savepath = strcat(pwd,'\Database_Reconstructions\Rec_',file_name,'.mat');
% % Save reconstruction
% save(output_savepath, 'output');
chunkData = zeros(258,506,Datasize-1);
NodeData = zeros(258,506,num_nodes);

for k = 1:num_nodes
    for kk = 1:Datasize-1
        chunkData(:,:,kk) = AllData{kk,k}.Statistical_Counts;
       
    end
    NodeData(:,:,k) = sum(chunkData,3);
    
    chunkData = zeros(258,506,Datasize-1);
end
output_savepath = strcat(pwd,'/Database_Reconstructions/Full_Rec_',outname,filename(1:end-5),'.mat');
save(output_savepath,'NodeData','enwind','AllFrame','X_rec','Y_rec','-v7.3');
disp('Saved!');

end


