function [ Frame ] = CalibrationDatasetLoad_Path( file_name, path)
%% This function loads calibration data contained in a .mat file 
%It takes as input the .mat file containing the acquired data (each row
%contains signal values from all the readout channels related to a single
%event)and the path of the folder containing it. Data are loaded into a 
%matrix called Frame, where values are saved with single precision

disp('Loading Dataset')
data = load(strcat(path,'Database\',file_name,'.mat'));
name_var = fieldnames(data);
expression = strcat('Frame = data.',name_var,';');
eval(expression{1});
Frame=single(Frame);
clc
end

