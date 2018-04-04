%%
clear all;
parameters.Run = ['079'];
parameters.Folder = 'C:\Users\Hanna\Documents\MATLAB\Research\2018-2\2-21\'; 
parameters.Sensitivity = 1;
parameters.Padding = 256;
if exist(strcat(parameters.Folder,'Figures')) == 0
    mkdir(strcat(parameters.Folder,'Figures'));
end

s3 = false;
shift = false;
IntN = 1;
conditions.G = 2.7; 
conditions.FixG = false; % replace fixG
conditions.FitPhase = false; % replace FitPhase
refshift = 1594.2;
underS = 2; 
shape = 2; % box: 0, decay: 1, echo: 2
data = csvread(strcat(parameters.Folder,'XCorrPP',parameters.Run,'.csv'),1,0);
%MergeTest3
%% Section into different pre-pulse scans
numInt = length(unique(data(:,19))); % Gives number of unique LC voltages/ pre pulse powers
IntTime = length(find(data(:,19) == data(1,19))); % Note that this works on the assumption that the number of steps is the same for each LC voltage
data = data(1+(IntN-1)*IntTime:IntN*IntTime,:); % This line is only meaningful if there is more than one LC voltage
%% Organize the data, add offsets, scale 
data = Sense(data); % Multiplies FWM and Linear signal by their respective sensitivities
tauSteps = length(unique(data(:,5)));
tSteps = length(unique(data(:,7)));
pos.t= unique(data(:,8));
pos.T = unique(data(:,4));
fprintf('hi')
[pos.tau,idxShort,idxLong] = unique(data(:,6)); % rename idxShort and idxLong to be more intuitive
offsetX = mean(data(:,9));
offsetY = mean(data(:,11));
procDataT = zeros(tSteps, 8, tauSteps);
procData = zeros(tSteps, 8, tauSteps);
for i = 1:tauSteps 
    tauRange = find(idxLong == i);
    tStart = find(pos.t == data(tauRange(1),8));
    tEnd = find(pos.t == data(tauRange(length(tauRange)),8));
    procDataT(tStart:tEnd,1,i) = data(tauRange,9)-offsetX;  % FWM X
    procDataT(tStart:tEnd,3,i) = data(tauRange,11)-offsetY; % FWM Y
    procDataT(tStart:tEnd,5,i) = data(tauRange,13); % Linear X
    procDataT(tStart:tEnd,7,i) = data(tauRange,15); % Linear Y
    procDataT(tStart:tEnd,2,i) = data(tauRange,10); % FWM X Sigma 
    procDataT(tStart:tEnd,4,i) = data(tauRange,12); % FWM Y Sigma
    procDataT(tStart:tEnd,6,i) = data(tauRange,14); % Linear X Sigma
    procDataT(tStart:tEnd,8,i) = data(tauRange,16); % Linear Y Sigma
end
%% Get phase shift from first (t = step 1, tau = step 1)
phase = zeros(2,1);
phase(2) = -atan2(procDataT(1,3,1),procDataT(1,1,1)); % FWM Signal phase
phase(1) = -atan2(procDataT(1,7,1),procDataT(1,5,1)); % Linear Signal phase
%%
    for I=1:tauSteps
        for J=1:tSteps
%     Real and imaginary parts of the expression R*e^(i phi)e^(i phi_0)
        procData(:,1,I) = procDataT(:,1,I)*cos(phase(2))-procDataT(:,3,I)*sin(phase(2)); % FWM X FWM Y % Real part
        procData(:,3,I) = procDataT(:,1,I)*sin(phase(2))+procDataT(:,3,I)*cos(phase(2));  %FWM X FWM Y % Imaginary part
        procData(:,5,I) = procDataT(:,5,I)*cos(phase(1))-procDataT(:,7,I)*sin(phase(1)); % Linear X Linear Y % Real part
        procData(:,7,I) = procDataT(:,5,I)*sin(phase(1))+procDataT(:,7,I)*cos(phase(1)); %Linear X Linear Y % Imaginary part
        end
    end
FWM.Real = procData(:,1,:); FWM.Imaginary = procData(:,3,:);
Linear.Real = procData(:,5,:); Linear.Imaginary = procData(:,7,:);
FWM.Complex = FWM.Real + j*FWM.Imaginary;
Linear.Complex = Linear.Real + j*Linear.Imaginary;
%% Make Time-Time Plots
timePlot(procData,parameters,pos,false)
%% Windowing that I do not understand
w1D = zeros(parameters.Padding,parameters.Padding);%tukeywin(2*parameters.Padding);
w1D(:,:) = 1;
w1D = tukeywin(2*tSteps); %length of tukeywin, r = 0.5
w1D = w1D(tSteps+1:end);
w1D2 = tukeywin(2*tauSteps); 
w1D2 = w1D2(tauSteps+1:end);
w2D = w1D*w1D2.';
w2D(:,:) = 1; 
w2Drep = repmat(w1D,1,tauSteps);
%{
% if (shift)
%     w2D(tSteps/2+1,tauSteps/2+1)=0;
% else
%     w2D(1,1)=0;
% end
    %}
%% Transform along t
FWM.tComplex= fft(squeeze(FWM.Complex),parameters.Padding,1);  
FWM.tReal = real(FWM.tComplex); FWM.tImaginary = imag(FWM.tComplex); % not phase shifted
Linear.tComplex = fft(squeeze(Linear.Complex),parameters.Padding,1); %w2D was taken out. w2D is needed if the photon echo is taken out. 
Linear.tReal = real(Linear.tComplex); Linear.tImaginary = imag(Linear.tComplex); % not phase shifted
%% Gets phase of data (tau = 1, t = all); 
phaseAll = -atan2(Linear.tImaginary(:,1), Linear.tReal(:,1)); % Uses only the first tau (because of photon echo shape) but would be ideal to use all the tau steps
phase2All = repmat(phaseAll,1,tauSteps);
FWM.tComplex = FWM.tReal.*cos(phase2All) - FWM.tImaginary.*sin(phase2All) + 1j*(FWM.tReal.*sin(phase2All) + FWM.tImaginary.*cos(phase2All));
%% Transform along tau
Linear.Transformed = fft(Linear.tComplex,parameters.Padding,2); 
FWM.Transformed = fft(FWM.tComplex,parameters.Padding,2);
%%
plot2DFWM = makeNice(FWM.Transformed,shift);
plot2DLinear = makeNice(Linear.Transformed,shift);
%% More things I do not understnad
    tStepSize = pos.t(2)-pos.t(1); % newVar
    tauStepSize = pos.tau(2)-pos.tau(1);
    tEnergy = 4.1357/tStepSize;
    tauEnergy = 4.1357/tauStepSize;
    BW = 4.1357/((pos.t(2)-pos.t(1))); % Bandwidth
    if shift
        freq = tEnergy*linspace(-0.5,0.5,parameters.Padding); 
        freq = freq+(freq(1)-freq(2))/2; % why the offset?
        freqtau = freq; % freq is really like freqt ..
    else % On currently. The vector has the same length as it would if shift = true, but is from 0 to 1 instead. 
            freq = BW*(linspace(1/(parameters.Padding+1),1,parameters.Padding));  %freq = tEnergy*(linspace(0,1,parameters.Padding)+1/parameters.Padding); % need to fix
            freqtau = -BW*(linspace(1,1/(parameters.Padding+1),parameters.Padding));
            plot2DFWM.Real = circshift(plot2DFWM.Real,1,1);
            plot2DFWM.Imaginary = circshift(plot2DFWM.Imaginary,1,1);
            plot2DFWM.Absolute = circshift(plot2DFWM.Absolute,1,1);
    end
    refshift = refshift+ underS*BW;
    freqref = freq+refshift; %takes the frequency and shifts it by 
    if s3
        freqtauref = -freqtau+2*refshift;
    else
        freqtauref = freqtau-refshift;
    end
parameters.Reference = freqref;
parameters.TauReference = freqtauref;
parameters.Frequency = freq;
%% Generate Absolute and Real Plots
generate2Dplot(plot2DFWM,parameters,false)
%% Making Appropriate Axes
[slices,projections,~,~] = sigWindow(parameters,plot2DFWM);  
%% Fit to Different Types of Lorentzians
[linewidths,fits] = fitLorentzian(slices,projections); 
%% Simulatenously fit the diagonal and anti-diagonal lineshapes for moderate inhomogeneity
[linewidths,fits] = fitSimultaneous(projections,linewidths,fits,slices,conditions,false); % I believe this is complete, but I have only tested it for fixG = false, fitPhase=false;
%% Plot Inhomogenous and Homogenous Slices with Simultaneous Fits
figure, 
subplot(1,2,1)
hold on
plot(projections.Homo+freq(21)+refshift,slices.Homo);
plot(projections.Homo+freq(21)+refshift,fits.SHomo);
set(gca,'FontSize',9);
xlabel('\omega_t (meV)')
axis square
title(strcat('\gamma = ',num2str(linewidths.Simultaneous),' meV'))
hold off
subplot(1,2,2)
hold on
plot(projections.Inhomo+freq(21)+refshift,slices.Inhomo);
plot(projections.Inhomo+freq(21)+refshift,fits.SInhomo);
set(gca,'FontSize',9);
xlabel('\omega_t (meV)')
axis square
hold off
%%



