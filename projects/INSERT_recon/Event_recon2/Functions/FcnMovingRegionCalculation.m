function [ index , data ] = FcnMovingRegionCalculation( coord, search_dim, Par  )

[ index_x ] = NN_finder(Par.x_knot,coord.x);
[ index_y ] = NN_finder(Par.y_knot,coord.y);

lim_x_inf = index_x - ceil(search_dim/2) + 1;
lim_x_sup = index_x + floor(search_dim/2);

lim_x_sup(lim_x_inf < 1) = search_dim;
lim_x_inf(lim_x_inf < 1) = 1;
lim_x_inf(lim_x_sup > length(Par.x_knot) ) = length(Par.x_knot) - search_dim +1;
lim_x_sup(lim_x_sup > length(Par.x_knot) ) = length(Par.x_knot);

lim_y_inf = index_y - ceil(search_dim/2) + 1;
lim_y_sup = index_y + floor(search_dim/2);


lim_y_sup(lim_y_inf < 1) = search_dim;
lim_y_inf(lim_y_inf < 1) = 1;
lim_y_inf(lim_y_sup > length(Par.y_knot) ) = length(Par.y_knot) - search_dim +1;
lim_y_sup(lim_y_sup > length(Par.y_knot) ) = length(Par.y_knot);

index.x = [lim_x_inf:lim_x_sup];
index.y = [lim_y_inf:lim_y_sup];

data.lim_x_inf = lim_x_inf;
data.lim_y_inf = lim_y_inf;
data.lim_x_sup = lim_x_sup;
data.lim_y_sup = lim_y_sup;
end