%% One for all script 

folder = uigetdir;
files = dir(fullfile(folder,'*.data'));
FilterSpec = '*.data';

for i = 1:length(files)

filename = files(i).name;
filepath = [files(i).folder,'\'];

[ NodeData, AllData ] = OneForAll(filepath,filename);

end