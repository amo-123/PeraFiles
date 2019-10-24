%% One for all script 

folder = uigetdir;
files = dir(fullfile(folder,'*.data'));
FilterSpec = '*.data';

for i = 1:1
filename = files(i).name;
filepath = [files(i).folder,'\'];

[ NodeData, AllData ] = OneForAll(filepath,filename,1,0);

end