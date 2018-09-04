function [fittedParameters,linewidths,fits,X,Ax] = fitSimultaneous(projections,linewidths,fits,slices,conditions,saveFig)
%[linewidths,fits] = fitSimultaneous(projections.Homo,projections.Inhomo,linewidths,fits)
% This works as of 3/14/2018 when fixG, fitPhase = false. Haven't tested
% for other conditions yet. 
peak = pkRatio(slices,linewidths,projections);
if conditions.FitPhase && conditions.FixG
     Ax = [real(slices.HomoC') real(slices.InhomoC') imag(slices.HomoC') imag(slices.InhomoC')];
     simultaneous = @(a,x) TO2X0hilbert(a,x,length(projections.Homo),length(projections.Inhomo),conditions.G);
     initialParameters = real([linewidths.Lorentzian peaks.Inhomo linewidths.GaussianAsyCenter 0 peaks.Homo 3 angle(sum(slices.InhomoC)) -0.58 2.55]); %linewidths.GaussianAsy]);
else
    Ax = [slices.Homo' slices.Inhomo'];
    X = [projections.Homo projections.Inhomo];
        if conditions.FixG
            initialParameters = real([linewidths.Lorentzian peaks.Inhomo linewidths.GaussianAsyCenter 0 peaks.Homo linewidths.LorentzianCenter]); %linewidths.GaussianAsy]);
            simultaneous = @(a,x) TO2X0fixG(a,x,length(projections.Homo),length(projections.Inhomo),conditions.G);
        else % turned on
            simultaneous = @(a,x) TO2X0(a,x,length(projections.Homo),length(projections.Inhomo),conditions.G);
            initialParameters = [linewidths.Lorentzian 2.8*peak.Inhomo linewidths.GaussianAsyCenter 0 peak.Homo linewidths.LorentzianCenter linewidths.GaussianAsy];
        end
end
opts.FunValCheck = 'off'; opts.MaxIter = 100;
[fittedParameters, ~, ~, covariance, ~] = nlinfit(X,Ax,simultaneous,initialParameters,opts);
linewidths.Simultaneous = fittedParameters(1);
linewidths.SimultaneousError = sqrt(diag(covariance));
fits.Simultaneous = simultaneous(fittedParameters,X); 
fits.SHomo = fits.Simultaneous(1:length(slices.Homo));
fits.SInhomo = fits.Simultaneous( (length(slices.Homo)+1):length(fits.Simultaneous));
fits.SimulCenter = fittedParameters(6);
%% Plot Inhomogenous and Homogenous Slices with Simultaneous Fits
% figure, 
% subplot(1,2,1)
% hold on
% plot(projections.Homo+freq(21)+refshift,slices.Homo);
% plot(projections.Homo+freq(21)+refshift,fits.SHomo);
% axis square
% hold off
% subplot(1,2,2)
% hold on
% plot(projections.Inhomo+freq(21)+refshift,slices.Inhomo);
% plot(projections.Inhomo+freq(21)+refshift,fits.SInhomo);
% axis square
% hold off
% if saveFig
%     print(strcat(parameters.Folder,'savedFigures\Run',num2str(parameters.Run),'LineshapeFits'),'-dpdf')        
% end
end

