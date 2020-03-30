%%Uniformity correction
function [K_interp] = UniformityCorrection(output,Par)
Counts = output.Statistical_Counts;

h1 = ones(5);
Centroid_Counts_filt =Par.sampling^2*conv2(Counts,h1,'same');

K = 1./Centroid_Counts_filt;

[X,Y]=meshgrid(linspace(min(Par.x_knot),max(Par.x_knot),size(Counts,2)),linspace(min(Par.y_knot),max(Par.y_knot),size(Counts,1)));
[Xq,Yq]=meshgrid(linspace(min(Par.x_knot),max(Par.x_knot),ceil(Par.cryst_lung_x/Par.pixel)),linspace(min(Par.y_knot),max(Par.y_knot),ceil(Par.cryst_lung_y/Par.pixel)));
K_temp = interp2(X,Y,K,Xq,Yq,'cubic');
K_interp = inpaint_nans(K_temp);

end