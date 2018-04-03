function [data] = Sense(data);
%
%%
sensitivity = [1 2 5 10 20 50 100 200 500];
% indices must be shifted since the first entry in sensitivity corresponds
% to the number 18. 
% I think this is missing some values. 
iXC = data(1,18)-10; 
iLinear = data(1,17)-10;
data(:,9:12) = data(:,9:12)*sensitivity(iXC);
%data(:,13:16) = data(:,9:12)*sensitivity(iLinear);
end

