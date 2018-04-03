function [linewidths,fits,func] = fitLorentzian(slices,projections)
% [linewidths,fits] = fitLorentzian(slices,projections.Homo,projections.Inhomo) 
% fitLorentzian fits the diagonal and cross diagonal slices to square root
% of a lorentzian, lorentzian, and asymmetric gaussian.

%% Square root of lorentzian
 clearvars fittedParameters covariance initialParameters
    func.LorentzianSquare = @(b,x) (b(2)*sqrt(b(1)./((2*pi)*((x-b(3)).^2+(b(1)).^2))));
    % b(1) = gamma, b(2) = integral of amplitude (i.e. simply the amplitude), b(3) = peak position
    initialParameters = [0.1,1,0]; 
    initialParameters(2) = -sum(slices.Homo)*(projections.Homo(2)-projections.Homo(1));
    [fittedParameters,~,~,covariance,~] = nlinfit(projections.Homo,transpose(slices.Homo),func.LorentzianSquare,initialParameters);
    linewidths.LorentzianSquare = fittedParameters(1);
    linewidths.LorentzianSquareError = sqrt(diag(covariance));
    fits.LorentzianSquare = func.LorentzianSquare(fittedParameters,projections.Homo);
%% 
 clearvars fittedParameters covariance initialParameters
    func.Lorentzian = @(b,x)(b(1)*b(2))./((2*pi)*((x-b(3)).^2+(b(1)).^2));
    % b(1) = HWHM, b(2) = integral of amp (i.e. simply the amp), b1(3) = peak pos
    %Fit homo line
    initialParameters= [0.1,1,0];
    initialParameters(2) = -sum(slices.Homo)*(projections.Homo(2)-projections.Homo(1));
    [fittedParameters,~,~,covariance,~] = nlinfit(projections.Homo,transpose(slices.Homo),func.Lorentzian,initialParameters);
    linewidths.Lorentzian = fittedParameters(1);
    linewidths.LorentzianCenter = fittedParameters(3);
    linewidths.LorentzianError = sqrt(diag(covariance));
    fits.Lorentzian = func.Lorentzian(fittedParameters,projections.Homo);
 %%
 clearvars fittedParameters covariance initialParameters
    func.GaussianAsymmetric = @(b,x)(b(2)*exp(-(x-b(3)).^2/(2*b(1)^2)).*(1+erf(b(4)*(x-b(3))/(b(1)))));
    lorentzianAsymmetric = @(b,x)((b(1)*b(2))./((2*pi)*((x-b(3)).^2+(b(1)).^2)).*(1+erf(b(4)*(x-b(3))/(sqrt(2)*b(1)))));
    %b(1) = HWHM, b(2) = integral of amp, b(3) = peak pos, b(4) = skew
    initialParameters = [1,1,0,0];
    initialParameters(2) = max(slices.Inhomo);
    [fittedParameters,~,~,covariance,~] = nlinfit(projections.Inhomo,transpose(slices.Inhomo),func.GaussianAsymmetric,initialParameters);
    linewidths.GaussianAsy = fittedParameters(1);
    linewidths.GaussianAsyCenter = fittedParameters(3);
    linewidths.GaussianAsyError = sqrt(diag(covariance));
    fits.GaussianAsy = func.GaussianAsymmetric(fittedParameters,projections.Inhomo);
  %%
 
end

