function [x_rec,y_rec,energy,error,Filt] = StatisticalMethod_for_LRF( Frame,Par,Tune,Filt )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
disp('Event Reconstruction')
%choose a sub dataset for Optical model 
Filt.subset_index = randsample(size(Frame,1),Tune.Num_rec);%random selection of a set number (Tune.Num_rec) of scintillation events. Their positions are reconstructed by using the LRFs (saved in Par.LRF).
Frame_sub = Frame(Filt.subset_index,:);

switch(Filt.Recon_Method)
    case 'ML'
        [x_rec,y_rec,energy,error] = ML_reconstruction( Par,Tune,Frame_sub );
    case 'LS'
        [x_rec,y_rec,energy,error] = LS_reconstruction( Par,Tune,Frame_sub );
    case 'WS'
        [Par.weights] = FcnWLSweightCalc(Frame,Par,Tune);
        [x_rec,y_rec,energy,error] = WSE_reconstruction( Par,Tune,Frame_sub );
end

clc
end

