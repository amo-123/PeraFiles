function [LRFs] = OpticalModelEngine( REC,Par,Tune, noise_on_data, fit_method)
%% This function estimates the Optical model (saved in the Light Response Function (LRF) map
%  A selected subset of events 
disp(['Optical Model Estimation'])
%% Initializations

%define subset values: 
%Frame is a matrix containing signals of the detection channels 
%related to the selected events only (randomly sampled and meeting the
%defined energy filters;
%x_rec and y_rec are the vectors of x and y coordinates of the selected
%events
Frame_rec=REC.Frame_rec;
x_rec=REC.x_rec;
y_rec=REC.y_rec;
biny = length(Par.y_knot);
binx = length(Par.x_knot);

%Light Response Function initialization
LRFs = cell(1,size(Frame_rec,2)); %cell array used to store the LRFs for all the photodetectors??
LRF=zeros(size(Frame_rec,2),biny,binx);% VARIABILE NON UTILIZZATA?

%Normalization of the dataset
Frame_rec_norm=Frame_rec./repmat(sum(Frame_rec,2),1,size(Frame_rec,2));
%sum(Frame_rec_2) provides a column vector whose elements represent the 
%overall signal for each gamma event; repmat copies this column vector a
%number of times equal to the number of photodetectors 

%Options depending on the fitting method: it allows to set the options for
%the parametric model to use for LRFs estimation (Gaussian/Bsplines)
[options]=OptionFitDistribution(fit_method);

%% Apply noise on data
if noise_on_data 
    %initialization of "exploration noise" to be applied on reconstructed
    %coordinates (noise values are randomly sampled from a normal
    %distribution with 0 mean and Par.sampling*Tune-dirty_noise standard
    %deviation
    noise_x=normrnd(0,Par.sampling*Tune.dirty_noise,[1 Tune.Num_rec]);
    noise_y=normrnd(0,Par.sampling*Tune.dirty_noise,[1 Tune.Num_rec]);
    
    %apply noise
    X_REC = x_rec+noise_x;
    X_REC(X_REC < 0) = 0;%if the "blurred" coordinates result to be negative, their value is set to zero (è il valore della coordinata più piccola possibile)
    X_REC(X_REC > Par.cryst_lung_x) = Par.cryst_lung_x;%if the blurred coordinates exceed the crystal length, their value is set to the crystal length (è il valore della coordinata più grande possibile)
    
    Y_REC = y_rec+noise_y;
    Y_REC(Y_REC < 0) = 0;
    Y_REC(Y_REC > Par.cryst_lung_y) = Par.cryst_lung_y;
else
    X_REC = x_rec;
    X_REC(X_REC < 0) = 0;
    X_REC(X_REC > Par.cryst_lung_x) = Par.cryst_lung_x;
    
    Y_REC = y_rec;
    Y_REC(Y_REC < 0) = 0;
    Y_REC(Y_REC > Par.cryst_lung_y) = Par.cryst_lung_y;
end

%% Estimate Optical Model through fitting for each detection channel
%  The i-th iteration cycle is used to estimate the LRF for the i-th
%  detection channel
for i=1:size(Frame_rec,2)
    F_REC=Frame_rec_norm(:,i);
    %the i-th photodetector center coordinates provide the position of the 
    %peak of the i-th LRF
    x0=Par.X_c(i);
    y0=Par.Y_c(i);

    if fit_method==1 %Gaussian Fitting
%         [fitresult] = GaussianFit2D(X_REC, Y_REC, F_REC, options, x0, y0);
        [fitresult] = GaussianFit2D_modified(X_REC, Y_REC, F_REC, options, x0, y0);
        LRFs{i}=fitresult;
    elseif fit_method==2
        n_div=7;
        LRFs{i} = Bspline_Fitting(Par,n_div,X_REC,Y_REC,F_REC,1,Tune);
    end    
end

%% Save Optical Model

clc
end

