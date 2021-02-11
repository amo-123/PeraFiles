%% One for all script 
function CallOneForAllScript(folder, uflag, Ufolder, EW, killmsk, outname)
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
else 
    Ufile = '~';
end

killdir = strcat(pwd,'/Database/',killmsk);

killchnl = load(killdir,'killch');

for i = 1:length(files)

filename = files(i).name;
filepath = [files(i).folder,'/'];

[ ~, ~ ] = OneForAll(filepath, filename, uflag, Ufolder, Ufile, EW, 1, killchnl.killch, outname);

end
toc;
end
