function [] = plot2D(signal,parameters,type,saveFig)
% Generate a 2D plot. The syntax is plot2D(signal,parameters,type,saveFig) 
% Signal - Fourier Transformed data
% Parameters - the parameters structure
% Title - self explanatory
% Savefig - True to save, false otherwise

figure();
    imagesc(parameters.Reference,parameters.TauReference,(signal)/parameters.Sensitivity);
    colormap('jet');
    if type == 1
        caxis([0,max(max(signal))]);
    end
    set(gca,'FontSize',14);
    set(gca,'YDir','normal')
    if type == 1
        title('Absolute')
    elseif type == 2
        title('Real')
    end
    line([parameters.Reference(1) parameters.Reference(parameters.Padding)], [-parameters.Reference(1) -parameters.Reference(parameters.Padding)],'Color','black');
    ylabel(strcat('\omega','_{\tau}',' (meV)'),'FontSize',16);
    xlabel('\omega_t (meV)','FontSize',16);
    colorbar;
if saveFig == true
    print(strcat(parameters.Folder,'\savedFigures\Run',num2str(parameters.Run)),'-fig')
end

end

