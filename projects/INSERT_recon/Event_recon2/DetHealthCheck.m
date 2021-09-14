function [ NodeFrame , NodeData, enwind, NF ] = DetHealthCheck(filepath,filename)

addpath(strcat(pwd,'/MATLAB_modified_Main_Split'));
addpath(strcat(pwd,'/MATLAB_modified_Main_Split/Functions'));
addpath(strcat(pwd,'/MATLAB_modified_Main_Split/Geometries'));
addpath(strcat(pwd,'/MATLAB_modified_Main_Split/Functions/dataFunctions'));


%% Initializations

%determine if the user wants to open only the last num_events
num_events=5000000; %number of events to display (comment this line to see all the events)
n_chan = 72;
%% Load            
%en_window_perc_width = 0.0865;
en_window_perc_width = 0.1;

%Datasize = 2;
enwind = zeros(20,3);
NF = zeros(20,72);

%f1 = figure();
%f2 = figure();
%f3 = figure('units','normalized','outerposition',[0 0 1 1]);
%f4 = figure();

ii = 1;
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
        NodeFrame = cell(num_nodes,1);
        NodeData = cell(num_nodes,1);
        
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
        
        
        rmpath(strcat(pwd,'/MATLAB_modified_Main_Split'));
        rmpath(strcat(pwd,'/MATLAB_modified_Main_Split/Functions'));
        rmpath(strcat(pwd,'/MATLAB_modified_Main_Split/Geometries'));
        rmpath(strcat(pwd,'/MATLAB_modified_Main_Split/Functions/dataFunctions'));
        
        
        
        for jj = 1:num_nodes
            %if jj == 14, continue, end
            goodimg = 0;
            while  goodimg == 0
            disp(strcat('Node #', num2str(jj),' Calculating...'));
            
            
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
            %addpath(strcat(pwd,'/Database'))
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
            
            %        addpath(strcat(pwd,'/Models_and_Corrections'));
            
            %        load('Optical_model_parameters.mat','Tune'); %gives 'Tune' structured variable as output
            
            % load dataset
            % Frame = CalibrationDatasetLoad(file_name);
            Frame = FRAME_NODE{jj};
        
            [~,ic1] = max( Frame, [], 2 );
            Frame = Frame( ic1>6, : );
           
        % Bad Event filter 
        if length(Frame) > 72
        B = sort(Frame,2);
        Diff = B(:,end) - B(:, end - 1) < mean(B(:,end) - B(:, end - 1))*2.8;
        Frame = Frame(Diff,:);
        
        M = mean(Frame(Frame>=30 & Frame<=500 )); % Kill channel 
        cut = 1.60;
        Mi = mean(Frame) >= M/cut & mean(Frame) <= M*cut;
        figure(1),  subplot(4,5,jj), plot(1:72,mean(Frame));
        hold on, plot(1:72,ones(1,72)*(M/cut),'LineWidth',2.5)
        hold on, plot(1:72,ones(1,72)*(M*cut),'LineWidth',2.5)
        Frame = Frame.*Mi;
        else
        end
       
            
            %figure(f1), title('unfiltered'), subplot(4,5,jj), imagesc(reshape(sum(Frame,1),[6,12])'), title('int2str(jj)'), axis image;
%             Frame_sub=Frame-600;
%             Frame_sub(Frame_sub<0)=0;
%             nx = 6;  ny = 12;
%             nv = size(Frame,1);
%             [~,ii_max] = max(Frame,[],2);
%             ix = mod( ii_max-1, nx ) + 1;
%             iy = fix( (ii_max-1) / nx ) + 1;
%             msk = zeros([nx,ny,nv]);
%             ix1=ix-2; ix1(ix1<1)=1; ix2=ix+2; ix2(ix2>nx)=nx;
%             iy1=iy-1; iy1(iy1<1)=1; iy2=iy+1; iy2(iy2>ny)=ny;
%             jx1=ix-1; jx1(jx1<1)=1; jx2=ix+1; jx2(jx2>nx)=nx;
%             jy1=iy-2; jy1(jy1<1)=1; jy2=iy+2; jy2(jy2>ny)=ny;
%             for iv=1:nv
%                 msk( ix1(iv):ix2(iv), iy1(iv):iy2(iv), iv )=1;
%                 msk( jx1(iv):jx2(iv), jy1(iv):jy2(iv), iv )=1;
%             end
%             msk = reshape( msk, [nx*ny,nv] )';
%             Frame = Frame_sub .* msk;
            %% ENERGY SPECTRUM
            % the EnergySpectrum function calculates the output energy (sum of the
            % signal from all channels) and plots the signal spectrum (uncalibrated energy spectrum
            %output.energy1 = EnergySpectrum(Frame, 0);
            
            %disp('%%%%%%%%%%%%%%%%%%%% PERA %%%%%%%%%%%%%%%%%%')
            %disp('>>> SELECTION of EVENTS to be reconstructed by statistical method <<<')
            %disp('Select energy values related to events to be reconstructed by the statistical method:')
            %disp('place the pointer on the lower limit of the range and click. Repeat the same for the upper limit.')
            
            %% Maual Peak
            
            %Energy for each of the acquired events is computed
            Energy = sum(Frame,2);
            %Energy histogram is computed (the number of bins is equal to the
            %sqrtof number of acquired events
            [y,x]=hist(Energy,floor(sqrt(size(Frame,1))));
            
            %The histogram is plotted: on the x-axis the energy amplitude is expressed
            %in ADC-channels
            if jj == 1, fe = figure('units','normalized','outerposition',[0 0 1 1]); end

            figure(fe);
            subplot(4,5,jj)
            hold on;
            plot(x,y,'LineWidth',2.5)
            axis([prctile(Energy,0.1) prctile(Energy,99.9) 0 1.1*max(y)])
            set(gca,'Fontsize',14,'Fontname','Arial','FontWeight','bold')
            xlabel('Signal [ADC channels]','Fontname','Arial','FontWeight','bold')
            ylabel('Counts [a.u]','Fontname','Arial','FontWeight','bold')
            title('Signal Spectrum','Fontname','Arial','FontWeight','bold')
            grid on
            
%             checkEW = 0;
%             while checkEW == 0
%             figure(fe);
%             [px,py]=ginput(1);
%             hold on;
%             plot(ones(size(1:py))*round((1-en_window_perc_width).*px),1:py,'LineWidth',2.5);
%             hold on;
%             plot(ones(size(1:py))*round((1+en_window_perc_width).*px),1:py,'LineWidth',2.5)
%              hold on;
%             plot(ones(size(1:py))*round(px),1:py,'LineWidth',2.5)
%             hold off; 
%             %checkEW = input('Is EW Correct Yes (1) or No (0)?');
%             %end
%             pause(1);
%             if jj == 16
%                 pause(0.5)
%             end
            [fitresult, ~] = fitRawFrame(x, y);
            FWHM  = 2*sqrt(log(2))*fitresult.c1;
            px = fitresult.b1;
            if fitresult.a1 > 10000
                fitresult.a1 = 5000;
            end
            %load('./EnergyWindows/EW_U01_L20191015.mat');
            %px = enwind(jj,3);
            figure(fe);
            subplot(4,5,jj)
            hold on
            plot(ones(size(1:fitresult.a1))*round(px),1:fitresult.a1,'LineWidth',1.5)
            hold on;
%             plot(px-(FWHM/2):px+(FWHM/2),ones(size(px-(FWHM/2):px+(FWHM/2)))*round(fitresult.a1/2),'LineWidth',1.5)
%             hold on;
            plot(ones(size(1:px))*round((1-en_window_perc_width).*px),1:px,'LineWidth',2.5);
            hold on;
            plot(ones(size(1:px))*round((1+en_window_perc_width).*px),1:px,'LineWidth',2.5)
            hold off; 
            xlim([0, 5e4]);
            pause(0.3);
            
            %k = 1; % select kth Channel to check single channel health
            
            %% Indi
            
            
            %ff = figure(); histogram(Frame(:,k));
            % Select the energy range for data energy windowing for image filtering)
            Filt.E_min=round((1-en_window_perc_width).*px);
            Filt.E_max=round((1+en_window_perc_width).*px);
            %Filt.E_min= min(px);%channels
            %Filt.E_max= max(px);%channels
            enwind(jj,1) = Filt.E_min;
            enwind(jj,2) = Filt.E_max;
            enwind(jj,3) = round(px);

            
            ch = sum(Frame);
            ch1 = ones(1,n_chan);
            if ( max(ch) > 0 )
                ii_chan = ( ch > 0 );
                ch1(ii_chan) = ch(ii_chan) / mean( ch(ii_chan) );
            end
            
            
            nrm = ones(size(Frame,1),1) * ch1;
            Frame = Frame ./ nrm;
            NF(jj,:) = ch1;
            
            %% CENTROID RECONSTRUCTION
            
            % CentroidReconstruction plots the image reconstructed with the
            % Centroid Method and returns the reconstructed coordinates of all the
            % events that follow the filter conditions
            %[output.x_rec_CM,output.y_rec_CM,output.Centroid_Counts,Filt.energy_window] = CentroidReconstruction(Frame,Par,Filt,output.energy1,0);
            %[output.x_rec_CM,output.y_rec_CM,output.Centroid_Counts,Filt.energy_window] = CentroidReconstructionKE(Frame,Par,Filt);
            disp('Centroid Recon')
            %[output.x_rec_CM,output.y_rec_CM,output.Centroid_Counts,Filt.energy_window] = CentroidRecon(Frame,Par,Filt);
            
            
            % Save initial dataset and work with a subset defined by the selected
            % energy window. Comment the two lines below to work with the wholw
            % dataset.
            %Frame_init = Frame;
            %Frame = Frame_init(Filt.energy_window,:);
            if jj == 1, f3 = figure('units','normalized','outerposition',[0 0 1 1]); end
            
             figure(f3), subplot(4,5,jj), imagesc(reshape(sum(Frame,1),[6,12])'), title(int2str(jj)), axis image; pause(0.9);
            
            %figure(f3), subplot(4,5,jj), imagesc([0 200], [0 400], output.Centroid_Counts'), title(int2str(jj)), axis image;
            
            %figure(f4), subplot(4,5,jj), bar(sum(Frame,1)), title(int2str(jj));
             %NodeFrame{jj} = Frame;
            
            %NodeData{jj} = output;
            pause(0.5);

%             goodimg = input('Is Recon Correct Yes (1) or No (0)?');
            goodimg = 1;
            %fn = [filename(1:end-5),'.png'];
            
%            close(fe);
            
            end
        
        end
        %saveas(f3,fn);

end


