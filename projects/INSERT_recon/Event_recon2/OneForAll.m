function [ NodeData, AllData ] = OneForAll(filepath,filename,EW,normflag,outname,doiflag)
% filepath: data file folder
% filename: data file name
% uflag: use flood spectra
% Ufilepath: flood file folder
% Ufilename: flood file name
% EW: Energy window file, if manual is required
% killchnl: matrix of pk targets
% outname: unique file identifier for output

    load(EW,'enwind');
    disp('Manual EW');


%% Main to Split Nodes
%Code to load .data files (INSERT files) and split the acquisitions from different nodes into FRAME_NODES variable

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

if doiflag == 1
    [Datasize] = MS_getfilesize(filename,filepath,num_events);
    AllData = cell(Datasize-1,20,4);
    disp('DOI Recon');
else
    [Datasize] = MS_getfilesize(filename,filepath,num_events);
    AllData = cell(Datasize-1,20);
end

FraLay = cell(1,4);

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
    NumFrame_NODE = FRAME_NODE(~cellfun('isempty',FRAME_NODE));
    num_nodes=length(NumFrame_NODE);
    
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

    for jj = n-num_nodes+1:n
            clear('output');
        disp(strcat('Node #', num2str(jj),' Calculating...'));

        %% LRFs for STATISTICAL RECONSTRUCTION
        % Set here the name of the file containing the LRFs (optical model) produced by 'LRF_estimation' script for the considered module.
        % The LRFs file is automatically saved by 'LRFs_estimation' in "LRFs"
        % folder.
        lrf_folder = './LRFs/';
        lrf_files = dir(fullfile(lrf_folder,'*.mat'));
        
        LRFs_filename = lrf_files(jj).name(1:end-4);
        
        
        %% BASELINE for MODIFIED CENTROID RECONSTRUCTION
        % Set the baseline value is subtracted to the signals of all the channels.
        % This is an arbitrary value and has to be defined by observing the result obtained in terms of reconstructed image employing the
        % Centroid Method: the baseline value should be set so that reconstructed events cover the whole FOV.
        % Typical values: from 400 to 600 ADC_channels.
        
        baseline = 400; % ADC_channels
        
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
%         for jj = 1:20
        Frame = FRAME_NODE{jj};

        %% ENERGY SPECTRUM
        % the EnergySpectrum function calculates the output energy (sum of the
        % signal from all channels) and plots the signal spectrum (uncalibrated energy spectrum
        %output.energy1 = EnergySpectrum(Frame, 0);
        
        %disp('%%%%%%%%%%%%%%%%%%%% PERA %%%%%%%%%%%%%%%%%%')
        %disp('>>> SELECTION of EVENTS to be reconstructed by statistical method <<<')
        %disp('Select energy values related to events to be reconstructed by the statistical method:')
        %disp('place the pointer on the lower limit of the range and click. Repeat the same for the upper limit.')
        
       
        Filt.E_min = enwind(jj,1);
        Filt.E_max = enwind(jj,2);
        
        
        %% Filter Channel Events
        [~,ic1] = max( Frame, [], 2 ); % Kill bright band phantoms 
        Frame = Frame( (ic1>6), : ); % Events in which the highest value is in channel 1 to 6 are removed.
       
        
        disp('Channel Normalisation');
        if exist('normflag','var') == 1
            load(normflag,'nrm');
            disp('Unifrom Normalisation');
            NF = ones(size(Frame,1),1) * nrm(jj,:);

            Frame = Frame ./ NF;
        end

        % Bad Event filter 
        if length(Frame) > 72
        B = sort(Frame,2);
        Diff = B(:,end) - B(:, end - 1) < mean(B(:,end) - B(:, end - 1))*2.8;
        Frame = Frame(Diff,:);
        
        M = mean(Frame(Frame>=30)); % Kill channel 
        cut = 1.6;
        Mi = mean(Frame) >= M/cut & mean(Frame) <= M*cut;
       
        Frame = Frame.*Mi;
        else
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
                if Mi(kk) == 0
                    LRFs{kk}.a = 0.02;
                    LRFs{kk}.b = 16;
                    LRFs{kk}.c = 16;
                    disp(strcat('LRF Node #', num2str(jj),' Channel #', num2str(kk), ' Killed.'));
                end
            end
        
            
            
            [Par.LRF] = LRF_Load(size(Frame,2),LRFs,Par);

           
            
            
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
            %[output.Statistical_Counts]=DisplayReconstruction(output,Par,Filt,Tune,0,0);
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
            
            %%  DOI calculation
            %filtraggio spaziale del Frame
            if doiflag == 1
                n_eventi=size(Frame,1);
                %number of groups
                n_groups=4;
                %% DOI Computation
                
                %here mean and std are evaluated starting from the LUT, then the log likelihood is computed and the Z groups are assigned where the likelihood is max
                
                %        Spotfiles = dir(fullfile('/SAN/inm/INSERT/amoINSERT/DOI/HSR/','*Frame*'));
                clear('SPOTS_fitting');
                Spotfiles = dir(fullfile('/media/ashley/My Passport/DOI/HSR/','*Frame*'));
                Spotfilename = strcat(Spotfiles(jj).folder,'/',Spotfiles(jj).name,'/SPOTs_fitting.mat');
                load(Spotfilename,'SPOTS_fitting');
                
                media=cell(n_groups,1);
                stdev=cell(n_groups,1);
                log_L=zeros(n_eventi, n_groups);
                Z_CLASS=zeros(n_eventi,1);
                
                for k=1:n_groups
                    disp(['Processing sigma and mean group ', num2str(k)])
                     parfor ch=1:n_chan
                        medias(:,ch)=feval(SPOTS_fitting{ch}.fitting_media{k},[output.x_rec', output.y_rec']);
                        stdevs(:,ch)= feval(SPOTS_fitting{ch}.fitting_stddev{k},[output.x_rec', output.y_rec']);
                    end
                    media{k} = medias;
                    stdev{k} = stdevs;
                end
                clear('medias');
                clear('stdevs');
                disp('DOI calculated');
                % use LUT to avoid this fitting ^
                
                for idx=1:n_eventi
                    for k=1:n_groups
                        log_L(idx,k)=-sum(ones(size(stdev{k}(idx,:)))*log(sqrt(2*pi))+log(stdev{k}(idx,:))+((Frame(idx,:)-media{k}(idx, :)).^2)./(2*stdev{k}(idx,:).*stdev{k}(idx,:)));
                    end
                end
                
                trovaM=@(vettore) find(vettore == max(vettore));
                for idx=1:n_eventi
                    Z_CLASS(idx)=trovaM(log_L(idx,:));
                end
                % Split Frame into 4 layers
                
                %%
                for k = 1:4
                    [MeshX,~] = meshgrid(Z_CLASS == k,1:72);
                    FraLay{k} = reshape(Frame(MeshX'),[],72);
                    Tune.Num_rec=size(FraLay{k},1);
                    %[output.x_rec,output.y_rec,output.energy,output.error,Filt] = StatisticalMethod(FraLay{k},Par,Tune,Filt );
                    
                    outputLay.energy = output.energy(Z_CLASS == k);
                    outputLay.x_rec_CM = output.x_rec_CM(Z_CLASS == k);
                    outputLay.y_rec_CM = output.y_rec_CM(Z_CLASS == k);
                    outputLay.x_rec = output.x_rec(Z_CLASS == k);
                    outputLay.y_rec = output.y_rec(Z_CLASS == k);
                    outputLay.error = output.error(Z_CLASS == k);
                    Par.pixel = 0.2;%mm
                    disp('ML DOI Recon');
                    [output.Statistical_Counts]=DisplayReconstruction(outputLay,Par,Filt,Tune,0,0);
                    
                    
                    %% SAVE RECONSTRUCTION
                    % Reconstructed event positions (X,Y), energy, reconstruction error are
                    % stored in the 'output' structured variable. This one is saved in
                    % Database_Reconstructions folder
                    AllData{ii,jj,k} = output;
                end
            else
                disp('ML Recon');
                [output.Statistical_Counts]=DisplayReconstruction(output,Par,Filt,Tune,0,0);
                AllData{ii,jj} = output;
            end
            
    end
        % <<<<<<< Updated upstream
        % =======
        rmpath(strcat(pwd,'/Geometries'));
        rmpath(strcat(pwd,'/Functions'));
        rmpath(strcat(pwd,'/Database'));
        
        addpath(strcat(pwd,'/MATLAB_modified_Main_Split'));
        addpath(strcat(pwd,'/MATLAB_modified_Main_Split/Functions'));
        addpath(strcat(pwd,'/MATLAB_modified_Main_Split/Geometries'));
        addpath(strcat(pwd,'/MATLAB_modified_Main_Split/Functions/dataFunctions'));
    
end

%% Save reconstruction

if doiflag == 1
    chunkData = zeros(258,506,Datasize-1,4);
    NodeData = cell(n,4);
    for k = 1:4
        for i = n-num_nodes+1:n
            for j = 1:Datasize-1
                
                if ~isempty(AllData{j,i,k}) == 1
                    chunkData(:,:,j,+k) = AllData{j,i,k}.Statistical_Counts;
                else
                    chunkData(:,:,j,kk) = chunkData(:,:,j,kk);
                end
            end
            NodeData{i,k} = sum(chunkData(:,:,:,k),3);
            
            chunkData = zeros(258,506,Datasize-1,4);
        end
        %NodeData = NodeData(~cellfun('isempty',NodeData));
    end
    
    output_savepath = strcat(pwd,'/Database_Reconstructions/DOI_Rec',outname,filename(1:end-5),'.mat');
    
    save(output_savepath,'NodeData','enwind','-v7.3');
    disp('Saved DOI!');
else
    chunkData = zeros(258,506,Datasize-1);
    NodeData = cell(n,1);
    
    for k = n-num_nodes+1:n
        for kk = 1:Datasize-1
            if ~isempty(AllData{kk,k}) == 1
                chunkData(:,:,kk) = AllData{kk,k}.Statistical_Counts;
            else
                chunkData(:,:,kk) = chunkData(:,:,kk);
            end
        end
        NodeData{k} = sum(chunkData,3);
        
        chunkData = zeros(258,506,Datasize-1);
    end
    
    %NodeData = NodeData(~cellfun('isempty',NodeData));
    
    output_savepath = strcat(pwd,'/Database_Reconstructions/Full_Rec_',outname,filename(1:end-5),'.mat');
    save(output_savepath,'NodeData','enwind','-v7.3');
    disp('Saved!');
end

end


