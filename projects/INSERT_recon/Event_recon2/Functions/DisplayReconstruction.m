function [Counts] = DisplayReconstruction( output,Par,Filt,Tune,isLRF,Unif_correction )
%% This function computes and displays the 2D histogram of reconstructed scintillation positions on the crystal XY plane
%  Reconstructed scintillation positions are filtered in order to select
%  only those events that follow filtering conditions related to energy,
%  reconstructed position and reconstruction error.Then, the 2D histogram 
%  of reconstructed scintillation positions on the XY plane is computed,
%  potentially corrected for uniformity and displayed.

% figure
    biny = ceil(Par.cryst_lung_y/Par.pixel);
    binx = ceil(Par.cryst_lung_x/Par.pixel);
    
    %The image displays only events that pass the following defined filters
    %related to energy, position reconstruction error and reconstructed 
    %position 
    
    F1=output.energy>Filt.E_min; 
    F2=output.energy<Filt.E_max;
    F3=(output.error>=Filt.error_min).*(output.error<=Filt.error_max);
    filt_x = (output.x_rec>0).*(output.x_rec<Par.cryst_lung_x);
    filt_y = (output.y_rec>0).*(output.y_rec<Par.cryst_lung_y);
    F4 = (filt_x.*filt_y);
    
    % Search for the indices of the events that match al the set filters
    Filt.kill_ML = find((F1.*F2.*F3.*F4)==1);
    
    % ??? Faccio diventare vuote le posizioni del vettore Filt.kill_ML in
    % corridpondenza delle quali sono salvati indici di eventi che superano
    % Tune.Num_rec?
    if isLRF == 1
        Filt.kill_ML(Filt.kill_ML>Tune.Num_rec)=[];
    end
    
    %Two "fake" events are added at the corners of the FOV to ensure the
    %histogram to get all the FOV (they are not considered in the output)
    y_image  = [output.y_rec(Filt.kill_ML) 0 Par.cryst_lung_y];
    x_image  = [output.x_rec(Filt.kill_ML) 0 Par.cryst_lung_x];

%    BY MICHELA: 
%     n_ev = size(x_image,2);
%     bin_tot = binx*biny;
%     sprintf('Numero di eventi per bin per avere una distribuzione uniforme: %d',n_ev/bin_tot)
     
    [Counts,Centers]=hist3([y_image' x_image'],[biny binx]);
    
    % The statistics of reconstructed events is shown
    stat.num_dati_tot = sum(sum(Counts));
    
    % The image borders are set to zero
    Counts(1,:)=0;Counts(:,1)=0;Counts(:,end)=0;Counts(end,:)=0;
    
    % If the 'K' field of Par is set and if 'Unif_correction' is equal to 
    % '1', the number of counts for each bin of the 2D histogram of
    % reconstructed position over the XY plane is corrected for uniformity.
    if any(strcmp('K',fieldnames(Par))) && Unif_correction
        %imagesc([0 Par.cryst_lung_x],[0 Par.cryst_lung_y],Counts.*Par.K)
        Counts = Counts.*Par.K;
    else
        %imagesc([0 Par.cryst_lung_x],[0 Par.cryst_lung_y],Counts)
    end
    
    %Settings for the 2D histogram of XY scintillation reconstructed
    %positions
%     set(gca,'Fontsize',15,'Fontname','Arial','FontWeight','bold')
%     xlabel('x [mm]')
%     ylabel('y [mm]')
%     hold on
%     colormap hot
%     axis equal
%     axis tight
%     colorb=colorbar;
%     set(colorb,'Fontsize',12);
%     set(colorb,'FontWeight','bold');
%     set(colorb,'Fontname','Arial');
%     drawnow;
end

