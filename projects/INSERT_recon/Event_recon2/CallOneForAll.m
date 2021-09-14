%% One for all script 
function CallOneForAll(folder, EW, nrm, outname,doiflag)
% folder: directory to get file from 
% uflag: use 1 for uniform flood spectra 
% Ufolder: Location of uni flood 
% EW: manual energy window file 
% outname: Filename Identifier 
tic;
% folder = uigetdir;
files = dir(fullfile(folder,'*.data'));
% FilterSpec = '*.data';

   % Ufiles = dir(fullfile(Ufolder,'*.data'));
  %  Ufile = Ufiles(1).name;
 %   Ufolder = [Ufolder,'/'];
    

%killdir = strcat(pwd,'/Database/',killmsk);

%killchnl = load(killdir,'killch');

for i = 6%:length(files)

filename = files(i).name;
filepath = [files(i).folder,'/'];

[ ~, ~ ] = OneForAll(filepath, filename, EW, nrm, outname,doiflag);

end
toc;
end
