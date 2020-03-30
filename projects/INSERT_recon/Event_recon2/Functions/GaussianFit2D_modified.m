function [fitresult, gof] = GaussianFit2D(X_REC, Y_REC, F_REC, opzioni, x0, y0)
%CREATEFIT(X_REC,Y_REC,F_REC)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : X_REC
%      Y Input : Y_REC
%      Z Output: F_REC
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 04-Dec-2014 14:44:02


%% Fit: 'untitled fit 1'.
X_REC =double(reshape(X_REC,length(X_REC),1));%this command allows to convert a row vector into a column one 
Y_REC =double(reshape(Y_REC,length(Y_REC),1));
F_REC =double(reshape(F_REC,length(F_REC),1));

[xData, yData, zData] = prepareSurfaceData( X_REC, Y_REC, F_REC );

% Set up fittype and options.
ft = fittype( 'a*exp(-(b*(x-x0)^2+c*(y-y0)^2)) + offset', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [opzioni.a_inf opzioni.b_inf opzioni.c_inf opzioni.offset_inf x0-6 y0-6];%options were set in OptionFitDistribution function
opts.StartPoint = [opzioni.a0 opzioni.b0 opzioni.c0 opzioni.offset0 x0 y0 ];
opts.Upper = [opzioni.a_sup opzioni.b_sup opzioni.c_sup opzioni.offset_sup x0+6 y0+6];

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );%[xData yData] = reconstructed positions of scintillation events;zData = signal values of a single readout channel (=> da un fotorivelatore) related to  all the events selected to perform LRF fitting
% 
% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, [xData, yData], zData );
% % h1 = plot( [xData, yData], zData );
% % hold on
% % plot( fitresult);
% % legend( h1, 'Gauss fit', 'Signal vs. X_R_E_C, Y_R_E_C', 'Location', 'NorthEast' );
% alpha(h,.5);
% % Label axes
% xlabel( 'X_R_E_C [mm]' );
% ylabel( 'Y_R_E_C [mm]' );
% zlabel( 'Signal [%]' );
% grid on
% view( 86.5, 16.0 );


