function [] = plotLinewidths(signal,linewidths,slices,fits,axes,parameters,saveFig)
% plotLinewidths calculates the horizontal shift off axis and plots
% lineshapes with fits
%% Get horizontal shift for homo
i = length(signal.Absolute)-axes.tauMax; 
diagonalPoint = [1+i,256-i]; 
maxPoint = [axes.tMax,axes.tauMax]; 
meV = parameters.Reference(2)-parameters.Reference(1);
roughShift = meV*(maxPoint(1) - diagonalPoint(1)); % gives detuning of the maximum from the diagonal line
fineShift = fits.SimulCenter; 
finalShift = roughShift + fineShift; 
%% Get right axes for inhomo
[itau,it] =ind2sub(size(signal.Absolute),find (signal.Absolute == max(slices.Inhomo)));
%%
figure, 
subplot(1,2,1)
    hold on
    plot(axes.Homo+roughShift,slices.Homo);
    plot(axes.Homo+roughShift,fits.SHomo);
    set(gca,'FontSize',9);
    xlabel('Emission Detuning from Diagonal (meV)') % should probably change this to something more like omega_t - omega diag???
    axis square
    title(strcat('\gamma = ',num2str(linewidths.Simultaneous),' meV .. shift = ',num2str(finalShift)))
    hold off
subplot(1,2,2)
    hold on
    plot(axes.Inhomo+parameters.Reference(it), slices.Inhomo);
    plot(axes.Inhomo+parameters.Reference(it), fits.SInhomo);
    set(gca,'FontSize',9);
    xlabel('\omega_t (meV)')
    axis square
    hold off
if saveFig
    print(strcat(parameters.Folder,'Figures\',num2str(parameters.Run),'_Lineshapes'),'-dpdf')        
end
end

