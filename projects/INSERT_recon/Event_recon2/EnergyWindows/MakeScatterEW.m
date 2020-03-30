%ScatterWindows


files = dir(fullfile(pwd,'*.mat'));

for i = 1:length(files)
    
files = dir(fullfile(pwd,'*.mat'));
filename = files(i).name;
filepath = [files(i).folder,'/'];


try 
load([filepath,filename],'enwind');
scatwind = zeros(20,2);
scatwind(:,2) = enwind(:,1);
diff = enwind(:,2) - enwind(:,1);
scatwind(:,1) = scatwind(:,2) - diff/2;
enwind = scatwind;
save([pwd,'\Scat',filename],'enwind')
catch 
    
end 

clear all

end