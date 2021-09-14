%% Make Energy windows

%folder = '/media/ashley/My Passport/Week_2/20190311/GeometricSensitivity/Det09/';
%folder = '/media/ashley/My Passport/MEASUREMENTS_London/20191001/Cylinder/';
%folder = '/media/ashley/My Passport/Week_1/20190306/Flood/';
folder = '/media/ashley/My Passport/London2021/Raw/20210908/Uni/';
%folder = '/media/ashley/My Passport/Week_2/20190313_part2/CylinderPhantom/';
%folder = '/media/ashley/My Passport/Week_1/20190308/forEnergyCalib/';
files = dir(fullfile(folder,'*.data'));

listwind = cell(length(files),1);
NF = cell(length(files),1);
Frame = cell(length(files),1);

for i = 1:length(files)

filename = files(i).name;
filepath = [files(i).folder,'/'];

%[ Frame{i} , ~, ~, ~ ] = DetHealthCheck(filepath,filename);

[ Frame{i} , ~, listwind{i}, NF{i} ] = DetHealthCheck(filepath,filename);
%close all

enwind = listwind{i};
save(['./EnergyWindows/EW_UF0',int2str(i), '_L20210908.mat'], 'enwind');
nrm = NF{i};
save(['./NrmFactors/NF_UF0',int2str(i), '_L20210908.mat'], 'nrm');

end
%%
% 
% enwind = listwind{2};
% save('./EnergyWindows/EW_U02_L20191216.mat', 'enwind');
% nrm = NF{2};
% save('./NrmFactors/NF_U02_L20191216.mat', 'nrm');
% % %%
% 
% enwind = zeros(20,2);
% for j = 17:20
%     
%   enwind(j,:) = listwind{j-1}(1,:);
% 
% end

%%
% PhotoPeaks = zeros(4,20);
% FWHM = zeros(4,20);
% Mu = zeros(4,20);
% fe = figure('units','normalized','outerposition',[0 0 1 1]);
% 
%   for i = 1:4
%       Data_Frame = Frame{i};
%       Data_EW = listwind{i};
%             figure(fe);
%             hold on;
%   for j = 1:20 
%             Energy = sum(Data_Frame{j},2);
%             px = Data_EW(j,3);
%             
%             %Energy histogram is computed (the number of bins is equal to the
%             %sqrtof number of acquired events
%             [y,x]=hist(Energy,floor(sqrt(size(Data_Frame{j},1))));
%             
%             %The histogram is plotted: on the x-axis the energy amplitude is expressed
%             %in ADC-channels
%             [fitresult, gof] = createFit(x, y);
%             FWHM(i,j) = 2*sqrt(log(2))*fitresult.c1;
%             Mu(i,j) = fitresult.b1;
%             
%             subplot(4,5,j);
%             hold on;
%             plot(x,y,'LineWidth',2.5)
%             axis([prctile(Energy,0.1) prctile(Energy,99.9) 0 1.1*max(y)])
%             set(gca,'Fontsize',14,'Fontname','Arial','FontWeight','bold')
%             xlabel('Signal [ADC channels]','Fontname','Arial','FontWeight','bold')
%             ylabel('Counts [a.u]','Fontname','Arial','FontWeight','bold')
%             title('Signal Spectrum','Fontname','Arial','FontWeight','bold')
%             grid on
%             hold on
%             plot(ones(size(1:fitresult.a1))*round(Mu(i,j)),1:fitresult.a1,'LineWidth',1.5)
%             hold on;
%             plot(Mu(i,j)-(FWHM(i,j)/2):Mu(i,j)+(FWHM(i,j)/2),ones(size(Mu(i,j)-(FWHM(i,j)/2):Mu(i,j)+(FWHM(i,j)/2)))*round(fitresult.a1/2),'LineWidth',1.5)
%             hold off; 
%             %checkEW = input('Is EW Correct Yes (1) or No (0)?');
%             %end
%             PhotoPeaks(i,j) = px;
%             [fitresult, gof] = createFit(x, y);
%             FWHM(i,j) = 2*sqrt(log(2))*fitresult.c1;
%             Mu(i,j) = fitresult.b1;
%             %pause(1);
%   end
%   end