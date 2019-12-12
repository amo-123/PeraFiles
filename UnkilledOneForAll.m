%% One for all script 
function UnkilledOneForAll(folder, uflag, Ufolder, EW, outname)
% folder: directory to get file from 
% uflag: use 1 for uniform flood spectra 
% Ufolder: Location of uni flood 
% EW: manual energy window file 
% outname: Filename Identifier 
tic;
% folder = uigetdir;
files = dir(fullfile(folder,'*.data'));
% FilterSpec = '*.data';

if uflag
    Ufiles = dir(fullfile(Ufolder,'*.data'));
    Ufile = Ufiles(1).name;
    Ufolder = [Ufolder,'/'];
end


for i = 1:length(files)

filename = files(i).name;
filepath = [files(i).folder,'/'];

[ NodeData, AllData ] = OneForAll(filepath, filename, 0, uflag, Ufolder, Ufile, EW, 1, outname);

end
toc;
end
