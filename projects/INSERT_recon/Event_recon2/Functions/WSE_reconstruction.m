function [ x_rec,y_rec,En_estimated , error] = WSE_reconstruction( Par,Filt,Frame)

%% Maximum Likelihood reconstruction
disp(['WSE reconstruction'])
% Filt.Num_rec=size(Frame,1);
x_rec=zeros(1,Filt.Num_rec);
y_rec=zeros(1,Filt.Num_rec);

LRF_temp=Par.LRF(:, ceil(Par.LRF_y_knot/Par.pixel_min) , ceil(Par.LRF_x_knot/Par.pixel_min) );
Par.sum_LRF_temp = reshape(sum(LRF_temp),size(LRF_temp,2),size(LRF_temp,3));
% Par.ln_LRF_temp = log(LRF_temp);

En_estimated=zeros(1,Filt.Num_rec);
Par.Search_dim_limit = 30; % dimensione massima della sottomatrice di ricerca nel metodo velocizzato
error = 1000*ones(1,Filt.Num_rec);

parfor i=2:Filt.Num_rec %parfor
%     flag = show_progress(i,round(Filt.Num_rec/100),Filt.Num_rec);
    % %%%%%%%%%%RICERCA CON GRIGLIA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%         [max_val,idx]=max(Frame(i,:));
%         coord = struct();
%         coord.x = Par.X_c(idx);
%         coord.y = Par.Y_c(idx);
%         search_dim.x = 1.25*Par.detect_dim + Par.dead_space;
%         search_dim.y = search_dim.x;
%     
%         search_dim_cost = 2;
%         end_condition = 0;
%         col_max=0;
%         row_max=0;
%     
%     
%         while ( end_condition ~= 1 )%(search_dim.x>=search_dim_cost*Par.pixel_min) || (search_dim.y>=search_dim_cost*Par.pixel_min) )
%     
%             [ search_dim ] = FncSearchDimRegulation( search_dim , Par );
%             [ grid_index , search_dim_new , end_condition, data ] = FcnRegionCalculation( coord, search_dim, Par  );
%     
%             % Metodo della massima verosimiglianza
%     
%             LRF_temp=Par.LRF(:, grid_index.y , grid_index.x );
%             %% implementazione Riccardo
%     
%             sum_LRF_temp = Par.sum_LRF_temp(grid_index.y , grid_index.x );
%             ln_LRF_temp = Par.ln_LRF_temp(:, grid_index.y , grid_index.x );
%             energ = sum(Frame(i,:));
%             N_estimated=repmat(energ,[size(LRF_temp,2) size(LRF_temp,3)])./sum_LRF_temp;
%             ln_L = sum(repmat(Frame(i,:)',[1 size(LRF_temp,2) size(LRF_temp,3)]).*ln_LRF_temp);
%             ln_L = energ*log(N_estimated)+reshape(ln_L,[size(LRF_temp,2) size(LRF_temp,3)]);
%             ln_L = ln_L - reshape(sum(LRF_temp,1),[size(LRF_temp,2) size(LRF_temp,3)]).*N_estimated;
%     
%             % Ricerca del Massimo
%             [max_val,idx]=max(ln_L(:));
%     
%             % %        caso in cui il minimo sia condiviso da più punti: bisogna
%             % %        distribuirlo
%             num_max=sum(sum(ln_L==max_val));
%             if num_max>1
%                 indici_cs = find(ln_L==max_val);
%                 pick_rand = ceil(num_max*rand(1));
%                 [row,col] = ind2sub(size(ln_L),indici_cs(pick_rand));
%             else
%                 [row,col] = ind2sub(size(ln_L),idx);
%             end
%     
%             if end_condition==1
%             data_norm = Frame(i,:)./sum(Frame(i,:));
%             data_model = LRF_temp(:,row,col)';
%             error(i) = sqrt(sum((data_norm-data_model).^2));
%             col_max=col;
%             row_max=row;
%             end
%       % update search window
%             coord.x = Par.LRF_x_knot(grid_index.x(col));
%             coord.y = Par.LRF_y_knot(grid_index.y(row));
%             search_dim.x = search_dim_new.x*search_dim_cost;
%             search_dim.y = search_dim_new.y*search_dim_cost;
%         end
%     
%         En_estimated(i) = N_estimated(row,col);
%         if(grid_index.x(col_max)<=1||grid_index.x(col_max)>=length(Par.LRF_x_knot)||grid_index.y(row_max)<=1||grid_index.y(row_max)>=length(Par.LRF_y_knot))%max sui bordi
%          x_rec(i) = coord.x;
%          y_rec(i) = coord.y;
%         else
%          grid_index.y=grid_index.y(row_max)+[-1,0,1];
%          grid_index.x=grid_index.x(col_max)+[-1,0,1];
%     
%          LRF_temp=Par.LRF(:, grid_index.y , grid_index.x );
%          sum_LRF_temp = Par.sum_LRF_temp(grid_index.y , grid_index.x );
%          ln_LRF_temp = Par.ln_LRF_temp(:, grid_index.y , grid_index.x );
%         % energ = sum(Frame(i,:));
%          N_estimated=repmat(energ,[size(LRF_temp,2) size(LRF_temp,3)])./sum_LRF_temp;
%          val_fit = sum(repmat(Frame(i,:)',[1 size(LRF_temp,2) size(LRF_temp,3)]).*ln_LRF_temp);
%          val_fit = energ*log(N_estimated)+reshape(val_fit,[size(LRF_temp,2) size(LRF_temp,3)]);
%          val_fit = val_fit - reshape(sum(LRF_temp,1),[size(LRF_temp,2) size(LRF_temp,3)]).*N_estimated;
%     
%          [delta_y,delta_x]=fitting_max_ML(val_fit);
%     
%     %      [delta_y,delta_x]=centroid_for_continuity(val_fit);
%          x_rec(i) = coord.x+delta_x*Par.pixel_min;
%          y_rec(i) = coord.y+delta_y*Par.pixel_min;
%     
%         end
    
    
    
    
    % %%%%%%%%%%RESEARCH WITH MOVING WINDOW%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [min_val,idx]=max(Frame(i,:));
    coord = struct();
    coord.x = Par.X_c(idx);
    coord.y = Par.Y_c(idx);
    search_dim = 9;
    
    end_condition = 0;
    col_max=0;
    row_max=0;
    
    
    while ( end_condition ~= 1 )%(search_dim.x>=search_dim_cost*Par.pixel_min) || (search_dim.y>=search_dim_cost*Par.pixel_min) )
        [ grid_index, ~ ] = FcnMovingRegionCalculation( coord, search_dim, Par  );

        % Metodo della massima verosimiglianza
        LRF_temp=Par.LRF(:, grid_index.y , grid_index.x );
        Weights_temp=Par.weights(:, grid_index.y , grid_index.x );
        %% implementazione Riccardo
        
        sum_LRF_temp = Par.sum_LRF_temp(grid_index.y , grid_index.x );
%         ln_LRF_temp = Par.ln_LRF_temp(:, grid_index.y , grid_index.x );
        energ = sum(Frame(i,:));
        N_estimated=repmat(energ,[size(LRF_temp,2) size(LRF_temp,3)])./sum_LRF_temp;
%         ln_L = sum(repmat(Frame(i,:)',[1 size(LRF_temp,2) size(LRF_temp,3)]).*ln_LRF_temp);
%         ln_L = energ*log(N_estimated)+reshape(ln_L,[size(LRF_temp,2) size(LRF_temp,3)]);
%         ln_L = ln_L - reshape(sum(LRF_temp,1),[size(LRF_temp,2) size(LRF_temp,3)]).*N_estimated;
        Chi = sum(  Weights_temp.*((repmat(Frame(i,:)',[1 size(LRF_temp,2) size(LRF_temp,3)])- permute(repmat(N_estimated,[1 1 36]),[3,1,2]).*LRF_temp).^2)  );
        Chi = permute(Chi,[2,3,1]);
        % Ricerca del Massimo
        [min_val,idx]=min(Chi(:));
        
        % %        caso in cui il minimo sia condiviso da più punti: bisogna
        % %        distribuirlo
        num_min=sum(sum(Chi==min_val));
        if num_min>1
            indici_cs = find(Chi==min_val);
            pick_rand = ceil(num_min*rand(1));
            [row,col] = ind2sub(size(Chi),indici_cs(pick_rand));
        else
            [row,col] = ind2sub(size(Chi),idx);
        end
        % if maximum at the border of the research window
        x_border=((col==1)||(col==size(Chi,2)));
        y_border=((row==1)||(row==size(Chi,1)));
        if x_border||y_border
            %if the maximum is at the edge of the crystal
            x_edges=((grid_index.x(col)==1)||(grid_index.x(col)==length(Par.LRF_x_knot)));
            y_edges=((grid_index.y(row)==1)||(grid_index.y(row)==length(Par.LRF_y_knot)));
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
        coord.x = Par.LRF_x_knot(grid_index.x(col));
        coord.y = Par.LRF_y_knot(grid_index.y(row));
        if end_condition
            data_norm = Frame(i,:)./sum(Frame(i,:));
            data_model = LRF_temp(:,row,col)';
            error(i) = sqrt(sum((data_norm-data_model).^2));
            col_max=col;
            row_max=row;
        end
    end
    
    En_estimated(i) = N_estimated(row,col);
    if(grid_index.x(col_max)<=1||grid_index.x(col_max)>=length(Par.LRF_x_knot)||grid_index.y(row_max)<=1||grid_index.y(row_max)>=length(Par.LRF_y_knot))%max sui bordi
        x_rec(i) = coord.x;
        y_rec(i) = coord.y;
    else
        y_cont=row_max+[-1,0,1];
        x_cont=col_max+[-1,0,1];
        val_fit = Chi(y_cont,x_cont);
        [delta_y,delta_x]=fitting_max_ML(val_fit);
        x_rec(i) = coord.x+delta_x*Par.pixel_min;
        y_rec(i) = coord.y+delta_y*Par.pixel_min;
    end
end



end

