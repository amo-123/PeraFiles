function [ Frame ] = CalibrationDatasetLoad( file_name)
%Summary of this function goes here
%   Detailed explanation goes here
disp('Loading Dataset')
data = load(strcat(pwd,'\Database\',file_name,'.mat'));
name_var = fieldnames(data);
expression = strcat('Frame = data.',name_var,';');
eval(expression{1});
Frame=single(Frame);
clc
end

