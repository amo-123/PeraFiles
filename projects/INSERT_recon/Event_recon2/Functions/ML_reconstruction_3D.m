function [ x_rec,y_rec,z_rec,En_estimated , error] = ML_reconstruction_3D( LUT_real,Par,Tune,Frame)

%% Maximum Likelihood reconstruction
disp(['ML reconstruction'])
% Tune.Num_rec=size(Frame,1);
x_rec=zeros(1,Tune.Num_rec);
y_rec=zeros(1,Tune.Num_rec);
z_rec=zeros(1,Tune.Num_rec);

limit_sampling = 0.2;

Par.LUT=LUT_real;
Par.sum_LUT = reshape(sum(Par.LUT,1),size(Par.LUT,2),size(Par.LUT,3),size(Par.LUT,4));
Par.ln_LUT = log(Par.LUT);

En_estimated=zeros(1,Tune.Num_rec);
Par.Search_dim_limit = 30; % dimensione massima della sottomatrice di ricerca nel metodo velocizzato
error = 1000*ones(1,Tune.Num_rec);

parfor i=1:Tune.Num_rec %parfor

    % %%%%%%%%%%RESEARCH WITH MOVING WINDOW%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [max_val,idx]=max(Frame(i,:));
    coord = struct();
    coord.x = Par.X_c(idx);
    coord.y = Par.Y_c(idx);
    search_dim = 9;
    
    end_condition = 0;
    col_max=0;
    row_max=0;
    
    
    while ( end_condition ~= 1 )%(search_dim.x>=search_dim_cost*Par.sampling) || (search_dim.y>=search_dim_cost*Par.sampling) )
        [ grid_index, ~ ] = FcnMovingRegionCalculation( coord, search_dim, Par  );

        %% Maximum Likelihood Method
        LUT_temp=Par.LUT(:, grid_index.y , grid_index.x,: );
        
        sum_LUT_temp = Par.sum_LUT(grid_index.y , grid_index.x ,:);
        ln_LUT_temp = Par.ln_LUT(:, grid_index.y , grid_index.x ,:);
        energ = sum(Frame(i,:));
        N_estimated=repmat(energ,[size(LUT_temp,2)  size(LUT_temp,3) size(LUT_temp,4)])./sum_LUT_temp;
        ln_L = sum(repmat(Frame(i,:)',[1 size(LUT_temp,2) size(LUT_temp,3) size(LUT_temp,4)]).*ln_LUT_temp);
        ln_L = energ*log(N_estimated)+reshape(ln_L,[size(LUT_temp,2) size(LUT_temp,3) size(LUT_temp,4)]);
        ln_L = ln_L - reshape(sum(LUT_temp,1),[size(LUT_temp,2) size(LUT_temp,3) size(LUT_temp,4)]).*N_estimated;
%         figure
%         subplot(2,2,1)
%         imagesc(reshape(ln_L(:,:,1),search_dim,search_dim))
%         colorbar
%         subplot(2,2,2)
%         imagesc(reshape(ln_L(:,:,2),search_dim,search_dim))
%         colorbar
%         subplot(2,2,3)
%         imagesc(reshape(ln_L(:,:,3),search_dim,search_dim))
%         colorbar
%         subplot(2,2,4)
%         imagesc(reshape(ln_L(:,:,4),search_dim,search_dim))
%         colorbar
%         drawnow
%         pause(1)
        % Ricerca del Massimo
        [max_val,idx]=max(ln_L(:));
        
        % %        caso in cui il minimo sia condiviso da più punti: bisogna
        % %        distribuirlo
        num_max=sum(sum(sum(ln_L==max_val)));
        if num_max>1
            indici_cs = find(ln_L==max_val);
            pick_rand = ceil(num_max*rand(1));
            [row,col,height] = ind2sub(size(ln_L),indici_cs(pick_rand));
        else
            [row,col,height] = ind2sub(size(ln_L),idx);
        end
        % if maximum at the border of the research window
        x_border=((col==1)||(col==size(ln_L,2)));
        y_border=((row==1)||(row==size(ln_L,1)));
        if x_border||y_border
            %if the maximum is at the edge of the crystal
            x_edges=((grid_index.x(col)==1)||(grid_index.x(col)==length(Par.x_knot)));
            y_edges=((grid_index.y(row)==1)||(grid_index.y(row)==length(Par.y_knot)));
            if x_edges||y_edges %the maximum is at (at least) one on the edges
                if x_edges && y_edges
                    %if the reconstruction coordinate is at the corner of the
                    %crystal --> stop the reconstruction
                    end_condition = 1;
                else
                    if not((x_edges && y_border)||(y_edges && x_border))
                    %if the coordinate is at one edge of the crystal and
                    %not at one border of the research area, stop the
                    %reconstruction
                    end_condition = 1;
                    end
                end                
            end
        else %the maximum is inside the research window --> stop condition
            end_condition = 1;
        end
        % update search coord
        coord.x = Par.x_knot(grid_index.x(col));
        coord.y = Par.y_knot(grid_index.y(row));
        coord.z = height;
        if end_condition
            data_norm = Frame(i,:)./sum(Frame(i,:));
            data_model = LUT_temp(:,row,col,height)';
            error(i) = sqrt(sum((data_norm-data_model).^2));
            col_max=col;
            row_max=row;
        end
    end
    
    En_estimated(i) = N_estimated(row,col);
    if(grid_index.x(col_max)<=1||grid_index.x(col_max)>=length(Par.x_knot)||grid_index.y(row_max)<=1||grid_index.y(row_max)>=length(Par.y_knot))%max sui bordi
        x_rec(i) = coord.x;
        y_rec(i) = coord.y;
        z_rec(i) = coord.z;
    else
        y_cont=row_max+[-1,0,1];
        x_cont=col_max+[-1,0,1];
        val_fit = ln_L(y_cont,x_cont,coord.z);
        [delta_y,delta_x]=fitting_max_ML(val_fit);
        x_rec(i) = coord.x+delta_x*Par.sampling;
        y_rec(i) = coord.y+delta_y*Par.sampling;
        z_rec(i) = coord.z;
        
%         y_cont=row_max+[-1,0,1];
%         x_cont=col_max+[-1,0,1];
%         val_fit = ln_L(y_cont,x_cont,:);
%        
%        sampling_z = 8/size(val_fit,3);%!!!!!!
%         Z0 = sampling_z/2:sampling_z:Par.cryst_lung_z;
%         
%          X0 = [-1,0,1]; Y0 = [-1,0,1];
%          X0=repmat(X0,size(val_fit,1),1,size(val_fit,3));
%          Y0=repmat(Y0',1,size(val_fit,2),size(val_fit,3));
%          Z0=repmat(reshape(Z0,1,1,length(Z0)),size(val_fit,1),size(val_fit,2),1);
%          limit_sampling = 0.2;
%         Xq = -1:limit_sampling:1; Yq = -1:limit_sampling:1; Zq = 0:0.5:Par.cryst_lung_z;
%         Xq = repmat(Xq,size(Yq,2),1,length(Zq));
%         Yq = repmat(Yq',1,size(Xq,2),length(Zq));
%         Zq = repmat(reshape(Zq,1,1,length(Zq)),size(Xq,1),size(Xq,2),1);
%         val_interp = interp3(X0,Y0,Z0,val_fit,Xq,Yq,Zq);
%         val_interp(:,:,end)=[];
%         val_interp(:,:,1)=[];
%         [max_val,idx]=max(val_interp(:));
%         [row,col,height] = ind2sub(size(val_interp),idx);
%         x_rec(i) = coord.x+Xq(1,col,1)*Par.sampling;
%         y_rec(i) = coord.y+Yq(row,1,1)*Par.sampling;
%         z_rec(i) =  Zq(1,1,height+1);
    end
end




end
