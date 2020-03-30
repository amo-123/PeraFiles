function [ LRF_par_spline ] = Bspline_Fitting(Par,n_div,X_REC,Y_REC,F_REC,immagine,Tune)

%Initialization 
x_rec=X_REC';
y_rec=Y_REC';
factor = 2;

[~,~,x_rec] = histcounts(x_rec,floor(Par.cryst_lung_x/Par.sampling/factor));
[~,~,y_rec] = histcounts(y_rec,floor(Par.cryst_lung_y/Par.sampling/factor));

F_REC=double(F_REC);
%sipm=1
%Frame =Frame./(repmat(sum(Frame,2),[1 size(Frame,2)]));%normalizzo
dim_x=size(Par.x_knot,2)/factor;
dim_y=size(Par.y_knot,2)/factor;

n_breaks_y = n_div;
n_breaks_x = n_div;
ky = 3;
kx = 3;
gain=0;%pixel%%%%%%%%%%%REMOVE
knotsx = augknt(linspace(1,dim_x,n_breaks_x+1),kx);
knotsy = augknt(linspace(1,dim_y,n_breaks_y+1),ky);


x_bin = 1:1:dim_x; y_bin = 1:1:dim_y;

%grid generation
grid=zeros(dim_y,dim_x);
cell_events=grid;%
for j=1:Tune.Num_rec%ingriglio i dati di Frame
  
    n_cell_events=cell_events(y_rec(j),x_rec(j));%num eventi cella fino a quel momento
    grid(y_rec(j),x_rec(j))=(F_REC(j)+grid(y_rec(j),x_rec(j))* n_cell_events)/(n_cell_events+1);%media
    cell_events(y_rec(j),x_rec(j))=n_cell_events+1;
end
grid=grid';

coef_spline=zeros(size(grid,1),length(knotsy)-ky);


prima_riga = 1;
ultima_riga = dim_x;
prima_colonna =1;
ultima_colonna = dim_y;

for i=1:size(grid,1)
    riga_griglia=grid(i,:);
    index_non_zeros=find(riga_griglia>0); 
    if(i>=prima_riga-gain&&i<=ultima_riga+gain&&length(index_non_zeros)>length(knotsy)-ky)%numero minimo di punti per soddisfare la regola di Schoenberg-Whitney 
         sp = spap2(knotsy,ky,y_bin(index_non_zeros),riga_griglia(index_non_zeros));
         coef_spline(i,:)=fnbrk(sp,'coefs');
    end
end
coefsy=coef_spline;
%    figure
%    imagesc(coefsy)

coefs=zeros(n_breaks_x+2,size(coefsy,2));
for j=1:size(coefsy,2)
    colonna_coef=coefsy(:,j);
    val_utili=find(colonna_coef>0);
    if (size(val_utili,1)>0)
        sp2 = spap2(knotsx,kx,val_utili,colonna_coef(val_utili)');
        coefs(:,j)= fnbrk(sp2,'coefs').';
    end
end
%  figure
%  imagesc(coefs) 
% % 
values = spcol(knotsx,kx,x_bin)*coefs*spcol(knotsy,ky,y_bin).'; 

% 
if immagine
figure
surf(values.','Linestyle','none');
drawnow
end
%  pause(2)
% close
% set(gca,'Fontsize',15,'Fontname','Arial','FontWeight','bold')
%     xlabel('x')
%     ylabel('y')
%     drawnow

LRF_par_spline.coefs=coefs;
LRF_par_spline.order_x=kx;
LRF_par_spline.order_y=ky;
LRF_par_spline.knotsx=knotsx;
LRF_par_spline.knotsy=knotsy;
LRF_par_spline.x_bin=x_bin;
LRF_par_spline.y_bin=y_bin;
end
