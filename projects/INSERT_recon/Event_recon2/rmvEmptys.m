% check empty 
clearvars;

folder = './Database_Reconstructions/Phantom/Milan/Sensitivity';
files = dir(fullfile(folder,'*.mat'));



for i = 1:length(files)
% [filename,filepath] = uigetfile([pwd,'\',FilterSpec], 'Select .data file', 'MultiSelect', 'off');

filename = files(i).name;
filepath = [files(i).folder,'\'];

data = open([filepath,filename]);
k = 1;
for j = 1: size(data.NodeData,3)
   if  sum(sum(data.NodeData(:,:,j))) ~= 0
       NodeData(:,:,k) = data.NodeData(:,:,j);
       k = k + 1;
       n = j;
   end

end
   fn = [filepath,filename]; 
   save(fn,'NodeData');
   clear NodeData;
end