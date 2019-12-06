%% One for all script 
function UnkillOneForAll(folder, uflag, Ufolder, EW, outname)
% folder: directory to get file from 
% uflag: use 1 for uniform flood spectra 
% Ufolder: Location of uni flood 
% EW: manual energy window file 
% outname: Filename Identifier 
tic;
% folder = uigetdir;
files = dir(fullfile(folder,'*.data'));
% FilterSpec = '*.data';

Ufiles = dir(fullfile(Ufolder,'*.data'));
Ufile = Ufiles(1).name

for i = 1:length(files)

filename = files(i).name;
filepath = [files(i).folder,'/'];

[ NodeData, AllData ] = OneForAll(filepath, filename, 0, uflag, Ufolder, Ufile, EW, outname);

end
toc;
end
