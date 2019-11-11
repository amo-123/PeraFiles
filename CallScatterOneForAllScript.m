%% One for all script 
function CallScatterOneForAllScript(folder)
% folder = uigetdir;
files = dir(fullfile(folder,'*.data'));
% FilterSpec = '*.data';

for i = 1:length(files)

filename = files(i).name;
filepath = [files(i).folder,'/'];

[ NodeData, AllData ] = ScatterOneForAll(filepath,filename,1,0);

end

end
