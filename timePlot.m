function [] = timePlot(procData,parameters,pos,FWM,saveFig)
%Generate time-time plots
%%
figure, hold on
subplot(1,2,1)
<<<<<<< HEAD
imagesc(pos.t,pos.tau,transpose(squeeze(procData(:,3,:))/parameters.Sensitivity))
ylabel('\tau (ps)'); xlabel('t (ps)');
axis square
subplot(1,2,2)
imagesc(pos.t,pos.tau,abs((FWM.Complex)/parameters.Sensitivity))
=======
imagesc(pos.t,pos.tau,transpose(squeeze(procData(:,3,:))/parameters.Sensitivity));
ylabel('\tau (ps)'); xlabel('t (ps)');
axis square
subplot(1,2,2)
imagesc(pos.t,pos.tau,abs(transpose(FWM.Complex)/parameters.Sensitivity));
>>>>>>> fbe82f3609e5727462a63207dc6ecc55ae0c5899
ylabel('\tau (ps)'); xlabel('t (ps)');
axis square
hold off
if saveFig
    print(strcat(parameters.Folder,'\Figures\',num2str(parameters.Run),'_TimeTime'),'-dpdf')
end
end

