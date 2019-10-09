function [x_rec,y_rec,energy,error,Filt] = StatisticalMethod( Frame,Par,Tune,Filt )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%disp'Event Reconstruction')
%choose a sub dataset for Optical model 
Frame_sub = Frame;

switch(Filt.Recon_Method)
    case 'ML'
        [x_rec,y_rec,energy,error] = ML_reconstruction( Par,Tune,Frame_sub );
    case 'LS'
        [x_rec,y_rec,energy,error] = LS_reconstruction( Par,Tune,Frame_sub );
    case 'WS'
        [Par.weights] = FcnWLSweightCalc(Frame,Par,Tune);
        [x_rec,y_rec,energy,error] = WSE_reconstruction( Par,Tune,Frame_sub );
    case 'M3'
        [x_rec,y_rec,energy,error] = ML_reconstruction_3D( Par,Tune,Frame_sub );
end

clc
end

