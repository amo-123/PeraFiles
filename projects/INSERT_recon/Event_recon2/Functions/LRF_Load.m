function [LRF] = LRF_Load(Num_channel,LRFs,Par)
%% This function samples previously loaded LRFs and saves them in LRF (3D matrix)
%   LRFs, loaded into a cell array (LRFs),are sampled over the XY crystal
%   plane according to Par.x_knot, Par.y_knot. The sampling result is saved
%   into LRF (3D matrix).

LRF=zeros(Num_channel,length(Par.y_knot),length(Par.x_knot)); %LRF = variable containing LRFs for all the readout channels (Num_channels). Each LRF is a 2D (x,y) map of values.

for i=1:Num_channel
    funzione=LRFs{i};
    LRF(i,:,:)=reshape(funzione([Par.x_coord',Par.y_coord']),length(Par.y_knot),length(Par.x_knot));%"funzione" viene valutata in corrispondenza delle coordinate Par.x_coord e Par.y_coord. Questi valori vengono riordinati in modo tale da crearsi una matrice 2D dei valori della LRF sul piano (x,y) del cristallo?  
end

end

