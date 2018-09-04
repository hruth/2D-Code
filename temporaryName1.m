%% Take a few diagonal slices detuned from the diagonal
A = [ 1:4; 5:8; 9:12; 13:16];
% one diagonal slice would be A(i,i) for all i
% another would be A(i,i+1) for all suitable i

%Main diagonal starts from 1,1

% diagSlice1 = zeros(10,length(plot2DFWM.Absolute));
diagSlice = zeros(10,length(plot2DFWM.Absolute));

% for i = 1:length(plot2DFWM.Absolute)
%     diagSlice1(1,i) = plot2DFWM.Absolute(i,i);
% end
% for i = 1:(length(plot2DFWM.Absolute)-1)
%     diagSlice1(2,i) = plot2DFWM.Absolute(i,i+1);
% end
idx1 = length(plot2DFWM.Absolute);
for i = 1:10;
idx = length(plot2DFWM.Absolute) - i + 1;
idx2 = sub2ind(size(plot2DFWM.Absolute),i:idx1, 1:idx);
diagSlice(i, 1:idx) = plot2DFWM.Absolute(idx2);
end

figure,hold on
for i = 1:10
    plot(diagSlice(i,:))
end
hold off

% can get (rough) off diagonal location from maxIndices

% figure,hold on
% plot(diagSlice1(1,:))
% plot(diagSlice1(2,:))
% hold off

%the line is defined by (1+i,256-i) for i in range 0 to 255
%%
figure,hold on
plot(simultaneous(fittedParameters,X))
plot(Ax)
hold off
%%
simultaneous = @(a,x) TO2X0(a,x,length(projections.Homo),length(projections.Inhomo),conditions.G);
Fit = simultaneous(fittedParameters,X);
figure,hold on
plot(projections.Homo,slices.Homo)
plot(projections.Homo,Fit(1:length(projections.Homo)))
hold off