function [slices,projections,tauMaxIndex,tMaxIndex] = sigWindow(parameters, nonlinearSignals)
%[slices,projections,tauMaxIndex,tMaxIndex] = sigWindow(parameters, nonlinearSignals) gets variables for creating a suitable window for
%the 2D plot as well as returning diagonal and cross diagonal slices. 
%%
 [M,I] = max(nonlinearSignals.Absolute(:)); [tauMaxIndex,tMaxIndex] = ind2sub(size(nonlinearSignals.Absolute),I); % Find peak location
homo = [1,1;parameters.Padding,parameters.Padding];
maxDiff = abs(tMaxIndex-tauMaxIndex);  maxSum = tMaxIndex+tauMaxIndex;
Length = parameters.Padding-maxDiff;
s = size(nonlinearSignals.Absolute);
slices.Homo = zeros(Length, 1); slices.HomoC = zeros(Length, 1);
% Xline = linspace(parameters.Reference(1),parameters.Reference(parameters.Padding),parameters.Padding);
% Yline = linspace(parameters.TauReference(1),parameters.TauReference(parameters.Padding),parameters.Padding);
% the coordinates [tau,t] will give you some pt on the two lines
if tMaxIndex <= tauMaxIndex %this one is false
    homo = [1+maxDiff 1; parameters.Padding parameters.Padding+maxDiff]   
    X=maxDiff; Y=0;
elseif tMaxIndex >= tauMaxIndex 
    homo = [1 1+maxDiff; parameters.Padding-maxDiff parameters.Padding]
    X=0; Y=maxDiff; % Gives the starting index of the cross-diagonal slice
end
  slices.Homo(1:Length) = nonlinearSignals.Absolute(sub2ind(size(nonlinearSignals.Absolute),X+(1:Length),Y+(1:Length)));
  slices.HomoC(1:Length) = nonlinearSignals.Real(sub2ind(size(nonlinearSignals.Absolute),X+(1:Length),Y+(1:Length))) + 1i*nonlinearSignals.Imaginary(sub2ind(size(nonlinearSignals.Absolute),X+(1:Length),Y+(1:Length)));
  line.Start = -abs(homo(1,1)-tauMaxIndex); line.End = abs(homo(2,1)-tauMaxIndex);
  projections.Homo = (parameters.Frequency(2)-parameters.Frequency(1))*linspace(line.Start,line.End,Length);
%   projections.Homo = (parameters.Frequency(2)-parameters.Frequency(1))*linspace(line.Start,line.End,Length);
%%
  %inhomo
    if maxSum >= parameters.Padding+1  
        inhomo = [parameters.Padding maxSum-parameters.Padding; maxSum-parameters.Padding parameters.Padding]
        Length = abs(maxSum-2*parameters.Padding)+1; X = parameters.Padding+1; Y=maxSum-parameters.Padding-1;
        slices.Inhomo = zeros(Length, 1); slices.InhomoC = zeros(Length, 1);
    else 
        inhomo = [maxSum-1 1; 1 maxSum-1];
        Length = maxSum-1; X = maxSum+1; Y = 0;
        slices.Inhomo = zeros(Length, 1); slices.InhomoC = zeros(Length, 1);
    end
        slices.Inhomo(1:Length) = nonlinearSignals.Absolute(sub2ind(s,X-(1:Length),Y+(1:Length)));
        slices.InhomoC(1:Length) = nonlinearSignals.Real(sub2ind(s,X-(1:Length),Y+(1:Length))) + 1i*nonlinearSignals.Imaginary(sub2ind(s,X-(1:Length),Y+(1:Length)));
        line.Start = -abs(inhomo(1,1)-tauMaxIndex); line.End = abs(homo(2,1)-tauMaxIndex);
        projections.Inhomo = (parameters.Frequency(2)-parameters.Frequency(1))*linspace(line.Start,line.End,Length);

end

