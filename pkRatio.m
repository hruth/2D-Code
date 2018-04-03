function [peak] = pkRatio(slices,linewidths,projections)
%Simple function to get pkRatios that are used as initial guesses in
%fitting functions.
%%
inputParameters = [linewidths.Lorentzian linewidths.GaussianAsy linewidths.LorentzianCenter 1];
homoLine = @(a,x) a(4)*abs((exp( ((abs(a(1))-1i*(x-a(3))).^2)/(2*abs(a(2))^2) ) .* (1+1i*erfi(1i*((abs(a(1)) - 1i*(x-a(3)))/(sqrt(2)*abs(a(2))))))) ./ (abs(a(2)).*(abs(a(1)) - 1i*(x-a(3))))); 
peak.Homo = max(max(slices.Homo))/max(max(homoLine(inputParameters,projections.Homo)));
peak.Inhomo = 2*max(max(slices.Inhomo));
end

