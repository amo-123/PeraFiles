function [fitresult, gof] = GaussFit(b, a,start)
%CREATEFIT(B,A)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : b
%      Y Output: a
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 25-Sep-2019 15:48:29


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( b, a );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
opts.StartPoint = [start 21116.26171875 1779.12820758099];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
% 
% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'a vs. b', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel b
% ylabel a
% grid on


