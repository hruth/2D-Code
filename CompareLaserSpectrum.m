%%
fileLocation = 'C:\Users\cundiff-lab\Documents\MATLAB\Hannas-Git-Repo\LaserSpectrum.csv';
% 2 to 2049
LaserSpectrum = csvread(fileLocation,0,0);
frequencyAxis = LaserSpectrum(:,1);
probe = LaserSpectrum(:,2);
pump1 = LaserSpectrum(:,4);
pump2 = LaserSpectrum(:,6);
figure,hold on
plot(frequencyAxis,probe)
plot(frequencyAxis,pump1)
plot(frequencyAxis,pump2)
hold off
%% PRESET 13 - 40 NM @ 740 
fileLocation = 'C:\Users\cundiff-lab\Documents\MATLAB\Hannas-Git-Repo\LaserSpectrum_Preset13.csv';
LaserSpectrum = csvread(fileLocation,0,0);
frequencyAxis = LaserSpectrum(:,1);
pulse = LaserSpectrum(:,2);
figure,hold on
title('Preset 13 - 40 nm @ 740')
plot(frequencyAxis,pulse)
xlabel('Wavelength (nm)')
ylabel('Intensity (arb. units)')
hold off 
%%
X = squeeze(procData(:,5,:));
Y = procData(:,7,:);
figure,hold on
subplot(2,1,1),plot(X(:,1))
subplot(2,1,2),plot(X(:,2))
hold off
%%
clearvars M
plot(X(:,1))
hfig = figure;
% f = getframe;
for k = 1:24
  plot(X(:,k))
  xlabel('t (s)')
%   legend(strcat('tau = ',num2str(2),'ps'))
  M(k) = getframe(hfig)
end
movie(gcf,M)
%%
Z = peaks(50);
hfig = figure;
hsurf = mesh(Z);
t=linspace(0,2*pi,16);
set(gca,'zlim',[-10 10]);
for i = 1:15,
      Zi = Z.*cos(t(i));
      set(hsurf,'zdata',Zi,'cdata',Zi);
      M(i) = getframe(hfig)
end;
figure
movie(gcf,M);
