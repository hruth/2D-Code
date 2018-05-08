function [] = timePlot(procData,parameters,pos,FWM,saveFig)
%Generate time-time plots
%%
figure, hold on
subplot(1,2,1)
contourf(pos.t,pos.tau,transpose(squeeze(procData(:,3,:))/parameters.Sensitivity),15,'edgecolor','none');
ylabel('\tau (ps)'); xlabel('t (ps)');
axis square
subplot(1,2,2)
imagesc(pos.t,pos.tau,abs(transpose(FWM.Complex)/parameters.Sensitivity),15,'edgecolor','none');
ylabel('\tau (ps)'); xlabel('t (ps)');
axis square
hold off
if saveFig
    print(strcat(parameters.Folder,'\savedFigures\TimeTimeRun',num2str(parameters.Run)),'-fig')
end
end

