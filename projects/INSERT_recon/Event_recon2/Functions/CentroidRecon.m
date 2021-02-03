function [x_rec,y_rec, Counts, energy_window,number_det_ch] = CentroidRecon( Frame, Par, Filt )
%Event reconstruction by centroid method

    %Define number of pixels in X and Y for image sampling
    biny = ceil( Par.cryst_lung_y / Par.pix(2) );
    binx = ceil( Par.cryst_lung_x / Par.pix(1) );
    
    % Definizione variabili
    E_min=Filt.E_min;              % energy filter values
    E_max=Filt.E_max;
    if isfield(Filt,'baseline')    % Baseline definition
        Baseline = Filt.baseline;
    else
        Baseline = mean( Frame(:) );  % channels
    end
    
    %Detection channels centers
    X_mat = repmat(Par.X_c,size(Frame,1),1);
    Y_mat = repmat(Par.Y_c,size(Frame,1),1);

    %Baseline subraction
    Frame_sub = Frame-Baseline;
    Frame_sub(Frame_sub<0)=0;

    %Neighbourhood mask (KE)
    if ( Filt.L_msk )
        Frame_sub = mask_frame( Frame_sub );
    end

    %Number of active detection channels per event, after baseline subtraction
    number_det_ch=sum(Frame_sub>0,2);

    %Centroid Formula
    sum_sub = sum(Frame_sub,2);
    x_rec = sum((Frame_sub.*X_mat),2) ./ sum_sub;
    y_rec = sum((Frame_sub.*Y_mat),2) ./ sum_sub;

    % Filters
    %energy window
    Energy = sum(Frame,2);
    F1 = Energy > E_min;
    F2 = Energy < E_max;
    
    %Coordinates reconstructed in the Field of View (FOV)
    x_crop = 7;
    F_tempX = ((x_rec>=x_crop)+(x_rec<=Par.cryst_lung_x))>1;
    F_tempY = ((y_rec>=0)+(y_rec<=Par.cryst_lung_y))>1;
    F_FOV = F_tempX.*F_tempY;
    
    %Filter dependent on the number of detection channels activated (influences only the image)
    F3 = number_det_ch <= size(Frame,2);
    F4 = number_det_ch >= 3;

    %Image filter (for Centroid reconstruction image)
    ii_filter = find( (F1.*F2.*F3.*F4.*F_FOV) == 1 );

    %new_filter_image = F1.*F2.*F3.*F4.*F_FOV;
    %filter_image = boolean(new_filter_image);

    y_image = y_rec(ii_filter,:)';
    x_image = x_rec(ii_filter,:)';
    
    %Two "fake" events added at the corners FOV to ensure the the whole FOV is included
    y_image  = [y_image 0 Par.cryst_lung_y];
    x_image  = [x_image 0 Par.cryst_lung_x];

    [Counts,Centers] = hist3([y_image' x_image'],[biny binx]);
    
    Counts(1,1) = Counts(1,1)-1;           % remove fake counts
    Counts(end,end) = Counts(end,end)-1;

    % Output reconstruction
    %Filter energy window (for LRF estimation)
    energy_window = find( (F1.*F2.*F_FOV) == 1 );
    y_rec = y_rec(energy_window,:)';
    x_rec = x_rec(energy_window,:)';

end
