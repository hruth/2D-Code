function [] = generate2Dplot(signal,parameters,saveFig)
% Generate a 2D plot. The syntax is plot2D(signal,parameters,type,saveFig) 
% Signal - Fourier Transformed data
% Parameters - the parameters structure
% Title - self explanatory
% Savefig - True to save, false otherwise
%%
figure,
% subplot(1,2,1)
imagesc(parameters.Reference,parameters.TauReference,(signal.Absolute)/parameters.Sensitivity);
colormap('jet');
caxis([0,max(max(signal.Absolute))]);
title('Absolute')
line([parameters.Reference(1) parameters.Reference(parameters.Padding)], [-parameters.Reference(1) -parameters.Reference(parameters.Padding)],'Color','black');
ylabel(strcat('\omega','_{\tau}',' (meV)')); xlabel('\omega_t (meV)');
set(gca,'YDir','normal');
colorbar; 
axis square
if saveFig
    print(strcat(parameters.Folder,'Figures\Run',num2str(parameters.Run),'Absolute'),'-dpdf')        
end

% subplot(1,2,2) 
figure,
imagesc(parameters.Reference,parameters.TauReference,(signal.Real)/parameters.Sensitivity);
colormap('jet');
title('Real')
line([parameters.Reference(1) parameters.Reference(parameters.Padding)], [-parameters.Reference(1) -parameters.Reference(parameters.Padding)],'Color','black');
ylabel(strcat('\omega','_{\tau}',' (meV)')); xlabel('\omega_t (meV)');
set(gca,'YDir','normal'); %'FontSize',14,'
colorbar;
axis square
if saveFig
    print(strcat(parameters.Folder,'Figures\Run',num2str(parameters.Run),'Real'),'-dpdf')        
end
end

