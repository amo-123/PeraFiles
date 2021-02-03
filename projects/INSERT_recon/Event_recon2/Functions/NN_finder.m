function [ index ] = NN_finder(array,coord)
coord = coord';
if length(coord) > 1
    index = sum((repmat(array,size(coord,1),1)-repmat(coord,1,size(array,2)))<0,2);
else
    index = sum((array-coord)<0);
end

% if index < length(array)
%     if int16(100*abs(array(index)-coord)/Par.pixel_min) ==  int16(100*abs(array(index+1)-coord)/Par.pixel_min)
%         if rand(1) < 0.5
%             index = index + 1;
%         end
%     end
% end

end

