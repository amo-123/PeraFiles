function [ x_rec,y_rec,En_estimated , error] = ML_reconstruction( Par,Filt,Frame)

%% Maximum Likelihood reconstruction
disp('ML reconstruction')
% Filt.Num_rec=size(Frame,1);
x_rec=zeros(1,Filt.Num_rec);
y_rec=zeros(1,Filt.Num_rec);

LRF_temp=Par.LRF(:, ceil(Par.y_knot/Par.sampling) , ceil(Par.x_knot/Par.sampling) );
Par.sum_LRF_temp = reshape(sum(LRF_temp),size(LRF_temp,2),size(LRF_temp,3));
Par.ln_LRF_temp = log(LRF_temp);

En_estimated=zeros(1,Filt.Num_rec);
Par.Search_dim_limit = 30; % dimensione massima della sottomatrice di ricerca nel metodo velocizzato
error = 1000*ones(1,Filt.Num_rec);

parfor i=1:Filt.Num_rec %parfor

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
        LRF_temp=Par.LRF(:, grid_index.y , grid_index.x );
        
        sum_LRF_temp = Par.sum_LRF_temp(grid_index.y , grid_index.x );
        ln_LRF_temp = Par.ln_LRF_temp(:, grid_index.y , grid_index.x );
        energ = sum(Frame(i,:));
        N_estimated=repmat(energ,[size(LRF_temp,2) size(LRF_temp,3)])./sum_LRF_temp;
        ln_L = sum(repmat(Frame(i,:)',[1 size(LRF_temp,2) size(LRF_temp,3)]).*ln_LRF_temp);
        ln_L = energ*log(N_estimated)+reshape(ln_L,[size(LRF_temp,2) size(LRF_temp,3)]);
        ln_L = ln_L - reshape(sum(LRF_temp,1),[size(LRF_temp,2) size(LRF_temp,3)]).*N_estimated;
        
        % Ricerca del Massimo
        [max_val,idx]=max(ln_L(:));
        
        % %        caso in cui il minimo sia condiviso da più punti: bisogna
        % %        distribuirlo
        num_max=sum(sum(ln_L==max_val));
        if num_max>1
            indici_cs = find(ln_L==max_val);
            pick_rand = ceil(num_max*rand(1));
            [row,col] = ind2sub(size(ln_L),indici_cs(pick_rand));
        else
            [row,col] = ind2sub(size(ln_L),idx);
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
        if end_condition
            data_norm = Frame(i,:)./sum(Frame(i,:));
            data_model = LRF_temp(:,row,col)';
            error(i) = sqrt(sum((data_norm-data_model).^2));
            col_max=col;
            row_max=row;
        end
    end
    
    En_estimated(i) = N_estimated(row,col);
    if(grid_index.x(col_max)<=1||grid_index.x(col_max)>=length(Par.x_knot)||grid_index.y(row_max)<=1||grid_index.y(row_max)>=length(Par.y_knot))%max sui bordi
        x_rec(i) = coord.x;
        y_rec(i) = coord.y;
    else
        y_cont=row_max+[-1,0,1];
        x_cont=col_max+[-1,0,1];
        val_fit = ln_L(y_cont,x_cont);
        [delta_y,delta_x]=fitting_max_ML(val_fit);
        x_rec(i) = coord.x+delta_x*Par.sampling;
        y_rec(i) = coord.y+delta_y*Par.sampling;
    end
end



end

