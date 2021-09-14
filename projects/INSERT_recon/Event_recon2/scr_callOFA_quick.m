
%folder = '/media/ashley/My Passport/Week_2/20190313_parct2/CylinderPhantom/';
folder = '/media/ashley/My Passport/London2021/Raw/20210818/LinX/';
%folder = '/media/ashley/My Passport/MEASUREMENTS_London/20191216/striatal/';
 uflag = 0;
 outname = 'tst_';
 killmsk = '20191216_killMask.mat';
Ufolder = '/media/ashley/My Passport/Week_1/20190306/Flood/';
 %EW = '/media/ashley/My Passport/TestLRF/PERA_PlanarReconstructionAlgorithm/PeraScripts/projects/INSERT_recon/Event_recon2/EnergyWindows/active/Universal_EW.mat';
%EW = [pwd,'/EnergyWindows/EW_U01_L20191216.mat'];
EW = [pwd,'/EnergyWindows/EW_UF02_L20210818.mat'];
doiflag = 0;
 nrm =  [pwd,'/NrmFactors/NF_UF01_L20210818.mat'];
CallOneForAll(folder, EW, nrm, outname,doiflag);
 
