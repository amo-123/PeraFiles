%% Make Energy windows

folder = '/media/ashley/My Passport/Week_2/20190311/GeometricSensitivity/';
%folder = '/media/ashley/My Passport/Week_2/20190312/GeometricCalibration/axial/';
%folder = '/media/ashley/My Passport/Week_1/20190306/Flood/';

files = dir(fullfile(folder,'*.data'));

listwind = cell(length(files),1);

for i = 1:length(files)

filename = files(i).name;
filepath = [files(i).folder,'/'];

[ ~ , ~, listwind{i} ] = DetHealthCheck(filepath,filename);
close all


end

enwind = listwind{1};
save('./EnergyWindows/EW_MGeoSens.mat', 'enwind');

%%

enwind = zeros(20,2);
for j = 1:20
    
  enwind(j,:) = listwind{j}(1,:);

end
  
  