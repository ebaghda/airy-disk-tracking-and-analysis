function [fitresult, gof] = FitMSD(dt, msd)
% This fits the MSDs for velocity
% X Input: dt
% Y Output: current_msd
% Output:
% fitresult : a fit object representing the fit.
% gof : structure with goodness-of fit info.
[xData, yData] = prepareCurveData(dt, msd);
%ft = fittype( '4*0.362*x+V*V*x*x', 'independent', 'x', 'dependent', 'y' ); %for water
ft = fittype( '4*0.182*x+V*V*x*x', 'independent', 'x', 'dependent', 'y' ); %for 20% IPA in water

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';
opts.StartPoint = 0.823457828327293;
[fitresult, gof] = fit( xData, yData, ft, opts );
end