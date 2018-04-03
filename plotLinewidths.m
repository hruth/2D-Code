Run = '013'
folder_h = 'D:\Hanna\Fits\2017\8-17-17\Mat Data\';
linewidthsFile = strcat(folder_h,'FittedLinewidths_Run',num2str(Run),'.mat');
load(linewidthsFile);

%%
numInt = 4;
figure
hold on
plot(sHomo(1:numInt),'*')
plot(aHomo(1:numInt),'*')
plot(abs(dHomo(1:numInt)),'*')
legend('Simultaneous Fit','Anti-Diagonal Fit','Diagonal Fit')
legend('Location','NorthOutside')
xlabel('Prepulse (uW)')
ylabel('linewidth (meV)')
title('Homogenous Linewidths')
axis square
hold off
%%
figure
hold on
plot(sInhomo(1:numInt),'*')
plot(aInhomo(1:numInt),'*')
plot(abs(dInhomo(1:numInt)),'*')
legend('Simultaneous Fit','Anti-Diagonal Fit','Diagonal Fit')
legend('Location','NorthOutside')
xlabel('Prepulse Iteration')
ylabel('linewidth (uW)')
title('Inhomogenous Linewidths')
axis square
hold off