clear all;

% File format:
% line 1: scan info (one tau step)
% lines 2-5: t sampled continuously at sampling rate (916 Hz for 0.5mm/s
% stage movement
%
% Generate array of data
% Find time 0 for each t scan
% Crop array accordingly
% Accountable: Eric Martin

% NEED TO FIX: work for 2 intensities

folder = 'D:\Eric\Data\2018\2018-09-02\';
% folder = 'D:\Eric\Data\2018\2018-08-31_cont\';
Run = '032';
%Run2 = '017';
Power = '1.0 uW'; %800 nW for dots
Tstep = 1;
imavg = 1; %number of images to average
sensitivity = 300; % intensities get large but are arbitrary, set this to get rid of big scaling in plots
%scaling = 10^5; % intensities get large but are arbitrary, set this to get rid of big scaling in plots
fits = true;
taustop = -1; %set -1 if not aborted
saveop = 0; % 0 for plot, 1 for save only, 2 for save and plot
hcn = 1.00029*1239842; %hc in air
cn = 1.00029*299792458;
refshift = 1586.3;%1534.91+0.63;%1595.7%4.2;%1628.9;%1627.3;

sampling = 916; %Hz from DAQ, should save this value in LV

%1696.5;%1644.9; %TMD monolayer
%1537.3 %1539.8;%1538.5; polariton (6/20, 5/23, 5/21)
%1626.7; dots
%1636.7; TMD
underS = 0;
taushift = 10.33;
s3 = false;
IntN = 1;
shape = 0; %box: 0, decay: 1, echo: 2
reflect = true;
plotextra = true;

% Set phase to some mean or 0 for determine phase later
%phase = [-1.4059;-2.7334]; %for 3/09: 005
%phase = [-1.9100;-1.6837]; %for 12/03: 054-8
phase = zeros(2,1); %set this by finding phase at zero

%S3 Run04
% maxamp = 648;%1584;
% maxamp2 = 2069;%1076;
% maxampr1 = 645;%1328;
% maxampi1 = 1891;%1053;
% maxampr2 = 596;%1580;
% maxampi2 = 2040;%978.5;
%S2+S1 (not right, check IntN = 4)
% maxamp = 1584;
% maxamp2 = 1076;
% maxampr1 = 1328;
% maxampi1 = 1053;
% maxampr2 = 1580;
% maxampi2 = 978.5;

%12/03: S3
% maxamp = 525.7;
% maxamp2 = 2737;
% maxampr1 = 496;
% maxampi1 = 484.7;
% maxampr2 = 2370;
% maxampi2 = 2664;

maxamp = 0;
maxamp2 = 0;
maxampr1 = 0;
maxampi1 = 0;
maxampr2 = 0;
maxampi2 = 0;

switch IntN
    case 1
        Power = ' 0 nW';
    case 2
        Power = ' 50 nW';
    case 3
        Power = ' 150 nW';
    case 4
        Power = ' 400 nW';
end

padi = 256;%1024; %number of points for fft (will pad with padi - scan length)
shift = false;

% probe_files =dir(strcat(folder,'Images',Run,'\Probe','*.bmp'));
% pump_files =dir(strcat(folder,'Images',Run,'\Pump','*.bmp'));
% numpos = length(probe_files)/imavg;
% thresh = 0.6;
% impr = zeros(480,744);
% impu = zeros(480,744);
% x = 1:744;
% y = 1:480;
% r = zeros(numpos,2); % r mean (1) and stddev (2)
% rtemp = zeros(imavg,1);
% rcut = 1.3/(6*0.004);
% rcut2 = 0.2/(6*0.004);
% 
% for I=1:numpos
%     for J=1:imavg
%         imprOrig = im2double(imread(strcat(folder,'Images',Run,'\',probe_files(imavg*(I-1)+J).name)));
%         impuOrig = im2double(imread(strcat(folder,'Images',Run,'\',pump_files(imavg*(I-1)+J).name)));
%         impr = squeeze(imprOrig(:,:,1));
%         impu = squeeze(impuOrig(:,:,1));
%         maxvalPr = max(max(impr));
%         maxvalPu = max(max(impu));
% 
%         imprNorm = impr./maxvalPr;
%         impuNorm = impu./maxvalPu;
% 
%         impr = impr.*fix(impr./maxvalPr + thresh);
%         impu = impu.*fix(impu./maxvalPu + thresh);
% 
%         Xpr = sum(impr,1);
%         Ypr = transpose(sum(impr,2));
%         PrSum = sum(Xpr);
%         Xpu = sum(impu,1);
%         Ypu = transpose(sum(impu,2));
%         PuSum = sum(Xpu);
% 
%         XPrMean = sum(Xpr.*x)/PrSum;
%         YPrMean = sum(Ypr.*y)/PrSum;
%         XPuMean = sum(Xpu.*x)/PuSum;
%         YPuMean = sum(Ypu.*y)/PuSum;
% 
%         rtemp(J) = sign(YPrMean-YPuMean).*sqrt((XPrMean-XPuMean)^2 + (YPrMean-YPuMean)^2);
%         plot(impr)
%         figure();
%         plot(impu)
%     end
%     r(I,1) = mean(rtemp(:));
%     r(I,2) = std(rtemp(:));
% end

data = csvread(strcat(folder,'XCorrPP',Run,'.csv'),1,0);
%data2 = csvread(strcat(folder,'XCorrPP',Run2,'.csv'),1,0);
%data(:,5) = roundn(data(:,5),-2);
%data(:,7) = roundn(data(:,7),-2);

%data = sortrows(data,[7 9]);

% determine # pump points, # LO points
% mirror positions = numpos = length(probe_files)
if s3
    temptau = 26;
else
    temptau = -1;
end
numtau = 0;
tempCapT = -50;
numCapT = 0;
TTime = 0;
stop = false;
widthEcho = 0; %really the half width
tempInt = 11;
numInt = 0;
IntTime = 0;
IntTimeCheck = 0; % error check for when shutter off had same LC value as 2nd
stgvel = data(1,7); % stage velocity
numt = length(data(1,:));

for J=1:length(data(:,1))/5
    I = (J-1)*5 + 1;
    %% Multiply data by sensitivity so values are in uV
    switch data(I,9)
        case 19
            sensitivity = 500;
        case 18
            sensitivity = 200;
        case 17
            sensitivity = 100;
        case 16
            sensitivity = 50;
        case 15
            sensitivity = 20;
        case 14
            sensitivity = 10;
        case 13
            sensitivity = 5;
        case 12
            sensitivity = 2;
        case 11
            sensitivity = 1;
    end
%     switch data(I,8)
%         case 19
%             data(I,13:16) = data(I,13:16)*500;
%         case 18
%             data(I,13:16) = data(I,13:16)*200;
%         case 17
%             data(I,13:16) = data(I,13:16)*100;
%         case 16
%             data(I,13:16) = data(I,13:16)*50;
%         case 15
%             data(I,13:16) = data(I,13:16)*20;
%         case 14
%             data(I,13:16) = data(I,13:16)*10;
%         case 13
%             data(I,13:16) = data(I,13:16)*5;
%         case 12
%             data(I,13:16) = data(I,13:16)*2;
%         case 11
%             data(I,13:16) = data(I,13:16)*1;
%     end
    %%
    % Determine data set to use for intensity iteration
    if (data(I,10) ~= tempInt)
        tempInt = data(I,10);
        numInt = numInt+1;
    end
    if numInt == 1
        IntTime = IntTime+5; %number of lines per iteration is 5
    end
    if numInt == 2
        IntTimeCheck = IntTimeCheck+5;
    end
    
    % Determine data set to use for T iteration
    if data(I,3) > tempCapT
        tempCapT = data(I,3);
        numCapT = numCapT+1;
    end
    if numCapT ==1
        TTime = TTime+5; %number of lines per iteration is 5
    end
    
    if s3
        if data(I,5) < temptau
            temptau = data(I,5);
            numtau = numtau+1;
        end
    else
        if data(I,5) > temptau
            temptau = data(I,5);
            numtau = numtau+1;
        end
    end
%     if shape == 2
%         if data(I,7) > tempt
%             tempt = data(I,7);
%             numt = numt+1;
%             if ~stop && numtau == 1
%                 widthEcho = widthEcho+1;
%             end
%         else
%             stop = true;
%         end
%     else
%         if data(I,7) > tempt && ~stop
%             tempt = data(I,7);
%             numt = numt+1;
%         else
%             stop = true;
%         end
%     end
end

%if IntTime > IntTimeCheck && numInt > 1
%    IntTime = IntTimeCheck
%end
%IntTime = IntTime/2;

if taustop ~= -1
    numtau = taustop;
end

%If first set int was set to same as next
%set data to only include point in IntN
data = data(1+(IntN-1)*IntTime:IntN*IntTime,:);

% for numpos
%startT = (Tstep-1)*tstep*numtau;
startT = (Tstep-1)*TTime;
Tpos = data(startT+1,4);

taupos = zeros(numtau,1);
for J = 1:numtau
    taupos(J,1) = 2*(data(1+(J-1)*5,5)-data(1,5))/0.1499;
end

procDataT = zeros(numt, 8, numtau);
procData = zeros(numt, 8, numtau);

Izero = zeros(numtau,1);

pWidth = 50; % # points on either side of zero to match pulse width
pcentWidth = 20;
pulse_range = zeros(pcentWidth*2+1,pWidth*2+1);

for I=1:numtau
    procDataT(:,1,I) = data(startT+(I-1)*5+2,:); %2X
    procDataT(:,3,I) = data(startT+(I-1)*5+3,:); %2Y
    procDataT(:,5,I) = data(startT+(I-1)*5+4,:); %1X
    procDataT(:,7,I) = data(startT+(I-1)*5+5,:); %1Y
    
    % Find time 0 using procDataT(1)
    % *TODO*: Find time 0 for all points after tau = 0 by fitting
    startPk = 1;%410;
    [M,Izero(I)] = max(procDataT(startPk:end,5,I));
    Izero(I) = Izero(I) + startPk;
    if s3
        Izero(I) = Izero(I)-round(taupos(I)*(cn*sampling)/(stgvel*0.002*10^12)); % minus becaus taupos is neg
    end
%     if (Izero(I) < 11200)
%         Izero(I) = 11350;
%     end
    
%     if I == 1
%         pulse_zero = procDataT(Izero(I)-pWidth:Izero(I)+pWidth,1,I);
%         pulse_zero = repmat(pulse_zero,pcentWidth*2+1,1);
%     else
%         for N = 1:pcentWidth*2+1
%             pulse_range(N,:) = Izero(I)-pWidth+N-pcentWidth-1:Izero(I)+pWidth+N-pcentWidth-1;
%         end
%         [Mn,Imin] = min(abs(sum(bsxfun(@minus,procDataT(transpose(pulse_range(:,:)),1,I),pulse_zero))));
%         Izero(I) = Izero(I) - pcentWidth + Imin + 1;
%     end

    if I == 1 && phase(1) == 0
        % TODO: fit data around zero with sinusoid
        if reflect
            phase(1) = -atan2(mean(procDataT(Izero(I)-5:Izero(I)+5,3,1)),mean(procDataT(Izero(I)-5:Izero(I)+5,1,1)));
        else
            phase(1) = -atan2(mean(procDataT(Izero(I)-5:Izero(I)+5,3,1)),mean(procDataT(Izero(I)-5:Izero(I)+5,1,1)));
        end
    end
    
    Izero(I) = 1;%Izero(I) - 300;
    
    procDataT(:,1,I) = circshift(procDataT(:,1,I),-Izero(I)+1); %shift zero to front
    procDataT(:,3,I) = circshift(procDataT(:,3,I),-Izero(I)+1); %shift zero to front
    procDataT(:,5,I) = circshift(procDataT(:,5,I),-Izero(I)+1); %shift zero to front
    procDataT(:,7,I) = circshift(procDataT(:,7,I),-Izero(I)+1); %shift zero to front
    
    %zero back of data (for decay window)
    procDataT(numt-2*Izero(I):numt,1,I) = 0;
    procDataT(numt-2*Izero(I):numt,3,I) = 0;
    procDataT(numt-2*Izero(I):numt,5,I) = 0;
    procDataT(numt-2*Izero(I):numt,7,I) = 0;
    
%     procData(:,1,I) = (procDataT(:,1,I)*phase(2,1)-procDataT(:,3,I)*phase(2,2))/norm(phase(2,:));
%     procData(:,2,I) = sqrt((procDataT(:,2,I)*phase(2,1)).^2+(procDataT(:,4,I)*phase(2,2)).^2)/norm(phase(2,:));
%     procData(:,3,I) = (procDataT(:,1,I)*phase(2,2)+procDataT(:,3,I)*phase(2,1))/norm(phase(2,:));
%     procData(:,4,I) = sqrt((procDataT(:,2,I)*phase(2,2)).^2+(procDataT(:,4,I)*phase(2,1)).^2)/norm(phase(2,:));
%     procData(:,5,I) = (procDataT(:,5,I)*phase(1,1)-procDataT(:,7,I)*phase(1,2))/norm(phase(1,:));
%     procData(:,6,I) = sqrt((procDataT(:,6,I)*phase(1,1)).^2+(procDataT(:,8,I)*phase(1,2)).^2)/norm(phase(1,:));
%     procData(:,7,I) = (procDataT(:,5,I)*phase(1,2)+procDataT(:,7,I)*phase(1,1))/norm(phase(1,:));
%     procData(:,8,I) = sqrt((procDataT(:,6,I)*phase(1,2)).^2+(procDataT(:,8,I)*phase(1,1)).^2)/norm(phase(1,:));

    %dlmwrite(strcat(folder,'XCorrPP_proc_',Run,'_',num2str(I),'.csv'),procData(:,:,I));
    
end
for I=1:numtau
%     for J=1:numt
%         if (J+I > numt)
%             procDataT(J,5,I) = 0;
%             procDataT(J,7,I) = 0;
%         end
%         if s3 && (J+I > numt)
%             procDataT(J,1,I) = 0;
%             procDataT(J,3,I) = 0;
%         end
%     end

    procData(:,5,I) = procDataT(:,1,I)*cos(phase(1))-procDataT(:,3,I)*sin(phase(1));
    %procData(:,6,I) = sqrt((procDataT(:,6,I)*phase(1,1)).^2+(procDataT(:,8,I)*phase(1,2)).^2)/norm(phase(1,:));
    procData(:,7,I) = procDataT(:,1,I)*sin(phase(1))+procDataT(:,3,I)*cos(phase(1));
    %procData(:,8,I) = sqrt((procDataT(:,6,I)*phase(1,2)).^2+(procDataT(:,8,I)*phase(1,1)).^2)/norm(phase(1,:));
    procData(:,1,I) = procDataT(:,5,I);
    procData(:,3,I) = procDataT(:,7,I);
end

% %cut off back of data
% numt = numt - max(Izero);
% procData = procData(1:numt,:,:);

tpos = linspace(0,numt*stgvel*0.002*10^12/(cn*sampling),numt);

%phase

%plots
% ept = steps;
% figure();
% if width(1,2) > 10
%     errorbar(taupos(2:ept,1),width(2:ept,1),width(2:ept,2),'ro-');
% else
%     errorbar(taupos(1:ept,1),width(1:ept,1),width(1:ept,2),'ro-');
% end
% title(strcat('Temporal response of spot size, P = ',Power,', LO Delay = ',num2str(Tpos,2)))
% ylabel('Sigma (um)')
% xlabel('Pump Delay (ps)')
% xlim([0,30])
% %ylim([620,670])
% set(gca,'box','off');
% ax1 = gca;
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none',...
%     'xticklabel',[],...
%     'ylim',[2.355*0.62,2.355*0.67]);
% ax2.YLabel.String = 'FWHM (um)';
% ax2.YLabel.Rotation = -90;

%%%

% for I = 1:numtau
%     procData(:,2,I) = procData(:,2,I)/max(procData(:,2,I));
% end

if (saveop == 0 || saveop == 2)
% ---TIME DOMAIN PLOTS---

% figure();
% contourf(tpos,taupos,transpose(squeeze(procData(:,5,:))/sensitivity),'edgecolor','none');
% title(strcat('Real part of SI, P = ',Power,', T Delay = ',num2str(Tpos,2)))
% ylabel('tau (ps)')
% xlabel('t (ps)')
% %xlim([0,6.3])
% %ylim([-6.3,0])
% 
% figure();
% contourf(tpos,taupos,transpose(squeeze(procData(:,5,:))/sensitivity),15,'edgecolor','none');
% title(strcat('Real part of SII, P = ',Power,', T Delay = ',num2str(Tpos,2)))
% ylabel('tau (ps)')
% xlabel('t (ps)')
% 
% figure();
% contourf(tpos,taupos,transpose(squeeze(procData(:,7,:))/sensitivity),'edgecolor','none');
% title(strcat('Imag part of SII, P = ',Power,', T Delay = ',num2str(Tpos,2)))
% ylabel('tau (ps)')
% xlabel('t (ps)')
% %xlim([0,6.3])
% %ylim([-6.3,0])

if plotextra
figure();
contourf(tpos,taupos,transpose(squeeze(abs(procData(:,5,:)+1j*procData(:,7,:)))/sensitivity),15,'edgecolor','none');
title(strcat('abs of SII, P = ',Power,', T Delay = ',num2str(Tpos,2)))
ylabel('tau (ps)')
xlabel('t (ps)')
% %xlim([0,6.3])
% %ylim([-6.3,0])
end

figure();
imagesc(tpos,taupos,transpose(squeeze(procData(:,1,:))/sensitivity));
set(gca,'YDir','normal')
%contourf(tpos,taupos,transpose(squeeze(procData(:,1,:))/sensitivity),15,'edgecolor','none');
colormap('redblue')
if s3
    title(strcat('Real part of S_{III}, P = ',Power,', T Delay = ',num2str(Tpos,2)))
else
    title(strcat('Real part of S_I, P = ',Power,', T Delay = ',num2str(Tpos,2)))
end
ylabel('tau (ps)')
xlabel('t (ps)')
%xlim([freqref(1),freqref(padi)]);

if plotextra

    
%xlim([freqref(1),freqref(padi)]);
%FINDME linear signal
figure();
imagesc(tpos,taupos,transpose(squeeze(abs(procData(:,1,:)+1j*procData(:,3,:)))/sensitivity));
set(gca,'YDir','normal')
%contourf(tpos,taupos,transpose(squeeze(abs(procData(:,1,:)+1j*procData(:,3,:)))/sensitivity),15,'edgecolor','none');
if s3
    title(strcat('Abs of S_{III}, P = ',Power,', T Delay = ',num2str(Tpos,2)))
    ylabel('T (ps)')
else
    title(strcat('Abs of S_I, P = ',Power,', T Delay = ',num2str(Tpos,2)))
    ylabel('\tau (ps)')
end
xlabel('t (ps)')
%xlim([freqref(1),freqref(padi)]);
end
end

%Window data
%tukey
w1D = kaiser(2*numt,25);%tukeywin(2*numt,1);
w1D = w1D(numt+1:end);
w1D2 = tukeywin(2*numtau);
w1D2 = w1D2(numtau+1:end);

% %exponential decay
% w1D = exp(-tpos/30);
% w1D2 = exp(taupos/15);
% 
w2D = w1D*w1D2.';
% w2D = zeros(numt,numtau);
% %no window
% w2D(:,:) = 1;

% if (shift)
%     w2D(numt/2+1,numtau/2+1)=0;
% else
%     w2D(1,1)=0;
% end

padit = 2^nextpow2(numt)*4;

sig2t = fft(w2D.*squeeze(procData(:,1,:)+1j*procData(:,3,:)),padit,1);
sig1t = fft(w2D.*squeeze(procData(:,5,:)+1j*procData(:,7,:)),padit,1);

addph = 0;
%phase2All = repmat(-atan2(imag(sig2t(:,1)),real(sig2t(:,1))),1,numtau);

%FINDME chirp correction
phase2All = -atan2(imag(sig2t),real(sig2t));
sig1t = (real(sig1t).*cos(phase2All+addph)-imag(sig1t).*sin(phase2All+addph))+1j*(real(sig1t).*sin(phase2All+addph)+imag(sig1t).*cos(phase2All+addph));

sig2 = fft(sig2t,padi,2);
sig1 = fft(sig1t,padi,2);

if (shift)
    reals2 = rot90(transpose(real(squeeze(fftshift(sig2)))),2);
    reals1 = rot90(transpose(real(squeeze(fftshift(sig1)))),2);
    imags2 = rot90(transpose(imag(squeeze(fftshift(sig2)))),2);
    imags1 = rot90(transpose(imag(squeeze(fftshift(sig1)))),2);
    abs2 = rot90(transpose(abs(squeeze(fftshift(sig2)))),2);
    abs1 = rot90(transpose(abs(squeeze(fftshift(sig1)))),2);
    %lin1 = rot90(transpose(abs(squeeze(fftshift(fft(fft(w2D.*squeeze(procData(:,1,:)+1j*procData(:,3,:)),padit,1),padi,2))))),2);
else
    reals2 = transpose(real(squeeze(sig2)));
    reals1 = transpose(real(squeeze(sig1)));
    imags2 = transpose(imag(squeeze(sig2)));
    imags1 = transpose(imag(squeeze(sig1)));
    abs2 = transpose(abs(squeeze(sig2)));
    abs1 = transpose(abs(squeeze(sig1)));    
%     reals2 = rot90(transpose(real(squeeze(sig2))),2);
%     reals1 = rot90(transpose(real(squeeze(sig1))),2);
%     imags2 = rot90(transpose(imag(squeeze(sig2))),2);
%     imags1 = rot90(transpose(imag(squeeze(sig1))),2);
%     abs2 = rot90(transpose(abs(squeeze(sig2))),2);
%     abs1 = rot90(transpose(abs(squeeze(sig1))),2);

%     reals1 = rot90(transpose(real(squeeze(fft(fft(w2D.*squeeze(procData(:,5,:)+1j*procData(:,7,:)),padit,1),padi,2)))),2);
%     imags1 = rot90(transpose(imag(squeeze(fft(fft(w2D.*squeeze(procData(:,5,:)+1j*procData(:,7,:)),padit,1),padi,2)))),2);
%     abs1 = rot90(transpose(abs(squeeze(fft(fft(w2D.*squeeze(procData(:,5,:)+1j*procData(:,7,:)),padit,1),padi,2)))),2);
%     lin1 = rot90(transpose(abs(squeeze(fft(fft(w2D.*squeeze(procData(:,1,:)+1j*procData(:,3,:)),padit,1),padi,2)))),2);
    %lin1 = rot90(transpose(abs(squeeze(fft(fft(w2D.*squeeze(procData(:,1,:)+1j*procData(:,3,:)),padit,1),padi,2)))),2);
end

% reals2 = w2D.*transpose(real(squeeze(fft(fft(procData(:,5,:)+1j*procData(:,7,:),padi,1),padi,3))));
% reals1 = w2D.*transpose(real(squeeze(fft(fft(procData(:,1,:)+1j*procData(:,3,:),padi,1),padi,3))));
% imags2 = w2D.*transpose(imag(squeeze(fft(fft(procData(:,5,:)+1j*procData(:,7,:),padi,1),padi,3))));
% imags1 = w2D.*transpose(imag(squeeze(fft(fft(procData(:,1,:)+1j*procData(:,3,:),padi,1),padi,3))));
% abs2 = w2D.*transpose(abs(squeeze(fft(fft(procData(:,5,:)+1j*procData(:,7,:),padi,1),padi,3))));
% abs1 = w2D.*transpose(abs(squeeze(fft(fft(procData(:,1,:)+1j*procData(:,3,:),padi,1),padi,3))));
freq = 4.13567/(tpos(2)-tpos(1))*linspace(0,1,padit)+1/padit;
if shift
    freqtau = -4.13567/((taupos(2)-taupos(1)))*linspace(-0.5,0.5,padi);
    freqtau = freqtau+(freqtau(1)-freqtau(2))/2;
else
    if s3
        freqtau = 4.13567/((taupos(2)-taupos(1)))*(linspace(0,1,padi)+1/padi);%+freq(1)-freq(2);
    else
        freqtau = -4.13567/((taupos(2)-taupos(1)))*(linspace(1,0,padi)+1/padi);
        % circshift SI tau axis so that DC is at end like all other axes
        reals1 = circshift(reals1,1,1);
        imags1 = circshift(imags1,1,1);
        abs1 = circshift(abs1,1,1);
        %lin1 = circshift(lin1,1,1);
    end
end

% if (saveop == 0 || saveop == 2)
% %-flipud(reals2)
% figure();
% contourf(freq,freqtau,(reals1)/sensitivity,15,'edgecolor','none');
% colormap('jet')
% title(strcat('Real part of SI, P = ',Power,', T Delay = ',num2str(Tpos,2),', run ',num2str(Run)))
% ylabel('wtau (meV)')
% xlabel('wt (meV)')
% axis equal;
% % xlim([-2,4]);
% % ylim([-3,3]);

if maxamp == 0
    maxamp = max(max(abs1));%23.97;38.49;
    %maxamp2 = max(max(abs2));
    maxampr1 = max(max((abs(reals1))));
    maxampi1 = max(max((abs(imags1))));
    %maxampr2 = max(max((abs(reals2))));
    %maxampi2 = max(max((abs(imags2))));
end

freqref = freq+refshift;
if s3
    refshift = 2*refshift- underS*(freqtau(padi)-freqtau(1)); %Include undersample
    freqtauref = -freqtau+refshift;
    %pltst = 1;
    %pltend = padit/10;
    [M,pltst] = min(abs(freqref(:)-freqtauref(1)/2));
    [M,pltend] = min(abs(freqref(:)-freqtauref(padi)/2));
else
    refshift = refshift+ underS*(freqtau(padi)-freqtau(1)); %Include undersample
    freqtauref = freqtau-refshift;
    %[M,pltst] = min(abs(freqref(:)+freqtauref(padi)));
    [M,pltend] = min(abs(freqref(:)+freqtauref(1)));
    pltst = 1;
end
pltend = pltend ;

spectfact = 1.4;
%%
figure();
imagesc(freqref(pltst:end),freqtauref,(abs1(:,pltst:end)/sensitivity));
%contourf(freqref,freqtauref,(abs1)/sensitivity,30,'edgecolor','none');%,15,'edgecolor','none');
colormap(jet);%viridis');
%caxis([0,(maxamp/sensitivity)]);%maxamp/5
set(gca,'FontSize',14);
set(gca,'YDir','normal')
if s3
    title(strcat('Absolute of S_{III}, T Delay = ',num2str(Tpos,2),' ps'))
    line([freqtauref(1) freqtauref(padi)], [2*freqtauref(1) 2*freqtauref(padi)],'Color','white','LineStyle',':');
    %Coupling line through 1645.09
    %line([freqref(1) freqref(padi)], [freqref(1)+1645.09 1645.09+freqref(padi)],'Color','white');
    %line([freqref(1) freqref(padi)], [freqref(1)+1644.42 1644.42+freqref(padi)],'Color','white');
    ylabel(strcat('\omega','_{T}',' (meV)'),'FontSize',16);
else
    title(strcat('Absolute of S_I, T Delay = ',num2str(Tpos,2),' ps'))%,', run ',num2str(Run)))
    line([-freqtauref(1) -freqtauref(padi)], [freqtauref(1) freqtauref(padi)],'Color','black','LineWidth',2,'LineStyle',':');
    ylabel(strcat('\omega','_{\tau}',' (meV)'),'FontSize',16);
end
xlabel('\omega_t (meV)','FontSize',16);
axis equal;
colorbar;
%%
xlim([freqref(pltst),freqref(end)]);
ylim([freqtauref(1),freqtauref(padi)]);

line([1540 1555], [-1545.84 -1545.84]*spectfact+1545.84*(spectfact-1),'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
line([1540 1555], [-1546.68 -1546.68]*spectfact+1545.84*(spectfact-1),'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
line([1540 1555], [-1547.39 -1547.39]*spectfact+1545.84*(spectfact-1),'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
line([1540 1555], [-1548.38 -1548.38]*spectfact+1545.84*(spectfact-1),'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
line([1540 1555], [-1549.08 -1549.08]*spectfact+1545.84*(spectfact-1),'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
line([1540 1555], [-1549.65 -1549.65]*spectfact+1545.84*(spectfact-1),'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');

line([1545.84 1545.84]*spectfact-1545.84*(spectfact-1), [-1540 -1555],'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
line([1546.68 1546.68]*spectfact-1545.84*(spectfact-1), [-1540 -1555],'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
line([1547.39 1547.39]*spectfact-1545.84*(spectfact-1), [-1540 -1555],'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
line([1548.38 1548.38]*spectfact-1545.84*(spectfact-1), [-1540 -1555],'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
line([1549.08 1549.08]*spectfact-1545.84*(spectfact-1), [-1540 -1555],'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
line([1549.65 1549.65]*spectfact-1545.84*(spectfact-1), [-1540 -1555],'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');

% line([1540 1555], [-1545.21 -1545.21]*spectfact+1545.21*(spectfact-1),'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
% line([1540 1555], [-1546.18 -1546.18]*spectfact+1545.21*(spectfact-1),'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
% line([1540 1555], [-1547.22 -1547.22]*spectfact+1545.21*(spectfact-1),'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
% line([1540 1555], [-1548.45 -1548.45]*spectfact+1545.21*(spectfact-1),'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
% line([1540 1555], [-1549.2 -1549.2]*spectfact+1545.21*(spectfact-1),'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
% line([1540 1555], [-1549.74 -1549.74]*spectfact+1545.21*(spectfact-1),'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
% 
% line([1545.21 1545.21]*spectfact-1545.21*(spectfact-1), [-1540 -1555],'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
% line([1546.18 1546.18]*spectfact-1545.21*(spectfact-1), [-1540 -1555],'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
% line([1547.22 1547.22]*spectfact-1545.21*(spectfact-1), [-1540 -1555],'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
% line([1548.45 1548.45]*spectfact-1545.21*(spectfact-1), [-1540 -1555],'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
% line([1549.2 1549.2]*spectfact-1545.21*(spectfact-1), [-1540 -1555],'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');
% line([1549.74 1549.74]*spectfact-1545.21*(spectfact-1), [-1540 -1555],'Color',[0,0,0]+0.3,'LineWidth',2,'LineStyle', '--');

% xlim([freqref(1),freqref(padit)]);
% ylim([-1643.5,-1638.5]);
% Export figure
%%
if s3
    saveas(gca,strcat(folder,num2str(Run),'_S3',num2str(IntN)),'epsc');
else
    saveas(gca,strcat(folder,num2str(Run),'_S1',num2str(IntN)),'epsc');
end
%export_fig(gca,strcat(folder,'S3',num2str(IntN),'2.eps'));

figure();
imagesc(freqref(pltst:pltend),freqtauref,(reals1(:,pltst:pltend))/sensitivity);
set(gca,'YDir','normal')
%contourf(freqref,freqtauref,(reals1)/sensitivity,30,'edgecolor','none');%,15,'edgecolor','none');
colormap('redblue');
caxis([-maxampr1/sensitivity,maxampr1/sensitivity])
set(gca,'FontSize',14);
if s3
    title(strcat('Real of S_{III}, P = ',Power,', T Delay = ',num2str(Tpos,2),' ps'))
    line([freqref(1) freqref(padi)], [freqtauref(1) freqtauref(padi)],'Color','white','LineStyle',':');
    ylabel(strcat('\omega','_{T}',' (meV)'),'FontSize',16);
else
    title(strcat('Real of S_I, P = ',Power,', T Delay = ',num2str(Tpos,2),' ps, run ',num2str(Run)))
    line([freqref(1) freqref(padi)], [freqtauref(padi) freqtauref(1)],'Color','white','LineStyle',':');
    ylabel(strcat('\omega','_{\tau}',' (meV)'),'FontSize',16);
end
xlabel('\omega_t (meV)','FontSize',16);
%axis equal;
colorbar;
xlim([freqref(pltst),freqref(pltend)]);
ylim([freqtauref(1),freqtauref(padi)]);
% xlim([1638.5,1643.5]);
% ylim([-1643.5,-1638.5]);
if s3
    saveas(gca,strcat(folder,num2str(Run),'_realS3',num2str(IntN)),'epsc');
else
    %saveas(gca,strcat(folder,num2str(Run),'_S1',num2str(IntN)),'epsc');
end

% % % if plotextra
% % % figure();
% % % imagesc(freqref,freqtauref,(imags1)/sensitivity);
% % % set(gca,'YDir','normal')
% % % %contourf(freqref,freqtauref,(imags1/max(max(imags1)))/sensitivity,30,'edgecolor','none');%,15,'edgecolor','none');
% % % colormap('redblue');
% % % caxis([-maxampi1/sensitivity,maxampi1/sensitivity])
% % % set(gca,'FontSize',14);
% % % if s3
% % %     title(strcat('Imag of S_{III}, P = ',Power,', T Delay = ',num2str(Tpos,2),' ps'))
% % %     line([freqref(1) freqref(padi)], [freqtauref(1) freqtauref(padi)],'Color','white','LineStyle',':');
% % %     ylabel(strcat('\omega','_{T}',' (meV)'),'FontSize',16);
% % % else
% % %     title(strcat('Imag of S_I, P = ',Power,', T Delay = ',num2str(Tpos,2),' ps, run ',num2str(Run)))
% % %     line([freqref(1) freqref(padi)], [freqtauref(padi) freqtauref(1)],'Color','white','LineStyle',':');
% % %     ylabel(strcat('\omega','_{\tau}',' (meV)'),'FontSize',16);
% % % end
% % % xlabel('\omega_t (meV)','FontSize',16);
% % % axis equal;
% % % colorbar;
% % % xlim([freqref(1),freqref(padit)]);
% % % % xlim([1638.5,1643.5]);
% % % % ylim([-1643.5,-1638.5]);
% % % end

% figure();
% imagesc(freqref,-freqtauref,flipud(reals2)/sensitivity);
% set(gca,'YDir','normal')
% %contourf(freqref,-freqtauref,flipud(reals2)/sensitivity,30,'edgecolor','none');%,15,'edgecolor','none');
% colormap('redblue');
% caxis([-maxampr2/sensitivity,maxampr2/sensitivity])
% set(gca,'FontSize',14);
% title(strcat('Real of S_{II}, P = ',Power,', T Delay = ',num2str(Tpos,2),' ps, run ',num2str(Run)))
% line([freqref(1) freqref(padi)], [-freqtauref(padi) -freqtauref(1)],'Color','white','LineStyle',':');
% ylabel(strcat('\omega','_{\tau}',' (meV)'),'FontSize',16);
% xlabel('\omega_t (meV)','FontSize',16);
% axis equal;
% colorbar;
% xlim([freqref(1),freqref(padi)]);

% if plotextra
% figure();
% imagesc(freq+refshift,-freqtauref,flipud(imags2/max(max(imags2)))/sensitivity);
% set(gca,'YDir','normal')
% %contourf(freq+refshift,-freqtauref,flipud(imags2/max(max(imags2)))/sensitivity,30,'edgecolor','none');%,15,'edgecolor','none');
% colormap('redblue');
% caxis([-maxampi2/sensitivity,maxampi2/sensitivity])
% set(gca,'FontSize',14);
% title(strcat('Imag of SII, P = ',Power,', T Delay = ',num2str(Tpos,2),', run ',num2str(Run)))
% line([freqref(1) freqref(padi)], [-freqtauref(padi) -freqtauref(1)],'Color','white','LineStyle',':');
% ylabel(strcat('\omega','_{\tau}',' (meV)'),'FontSize',16);
% xlabel('\omega_t (meV)','FontSize',16);
% axis equal;
% colorbar;
% xlim([freqref(1),freqref(padi)]);
% end

% if ~s3
% figure();
% imagesc(freq+refshift,-freqtauref,flipud(abs2/sensitivity).^2);
% set(gca,'YDir','normal')
% %contourf(freq+refshift,-freqtauref,flipud(abs2/sensitivity).^2,30,'edgecolor','none');
% colormap('inferno');
% caxis([0,(maxamp2/sensitivity)^2]);%maxamp2/5
% set(gca,'FontSize',14);
% title(strcat('Absolute part of S_{II}, P = ',Power,', T Delay = ',num2str(Tpos,2),' ps'))
% line([freqref(1) freqref(padi)], [-freqtauref(padi) -freqtauref(1)],'Color','white','LineStyle',':');
% ylabel(strcat('\omega','_{\tau}',' (meV)'),'FontSize',16);
% xlabel('\omega_t (meV)','FontSize',16);
% axis equal;
% colorbar;
% xlim([freqref(1),freqref(padi)]);
% % xlim([-2,4]);
% % ylim([-3,3]);
% saveas(gca,strcat(folder,num2str(Run),'_S2',num2str(IntN)),'epsc');
% end

% if ~s3
% figure();
% imagesc(freq+refshift,-freqtauref,flipud(reals2+flipud(reals1))/sensitivity);
% set(gca,'YDir','normal')
% %contourf(freq+refshift,-freqtauref,flipud(reals2+flipud(reals1))/sensitivity,30,'edgecolor','none');
% colormap('viridis');
% set(gca,'FontSize',14);
% title(strcat('Correlation Spectrum, P = ',Power,', T Delay = ',num2str(Tpos,2),' ps'))
% %line([freqref(1) freqref(padi)], [-freqtauref(padi) -freqtauref(1)],'Color','black');
% ylabel('w_{\tau} (meV)')
% xlabel('w_t (meV)')
% axis equal;
% colorbar;
% xlim([freqref(1),freqref(padi)]);
% end

% Plot phase
figure();
surf(tpos,taupos,transpose(squeeze(atan2(procData(:,7,:),procData(:,5,:)))),'edgecolor','none');
title(strcat('SI Phase, P = ',Power,', T Delay = ',num2str(Tpos,2),' ps'))
ylabel('tau (ps)')
xlabel('t (um)')

% figure();
% contourf(freq+refshift,freqtau-refshift,atan2(imags1,reals1)/pi,10,'edgecolor','none');
% %title(strcat('SI Phase, P = ',Power,', T Delay = ',num2str(Tpos,2)))
% colormap('jet');
% set(gca,'FontSize',14);
% ylabel(strcat('\omega','_{\tau}',' (meV)'),'FontSize',16);
% xlabel('\omega_t (meV)','FontSize',16);
% %axis equal;
% xlim([1642.5,1644]);
% ylim([-1644.8,-1643.3]);
% caxis([-1,1]);

% figure();
% contourf(freq+refshift,freq-refshift,(atan2(imags1,reals1).*1./(1+exp(-25*(abs1/maxamp-0.2))))/sensitivity,20,'edgecolor','none');%,15,'edgecolor','none');
% colormap('redblue');
% set(gca,'FontSize',14);
% %title(strcat('Absolute of SI, P = ',Power,', T Delay = ',num2str(Tpos,2),', run ',num2str(Run)))
% ylabel(strcat('\omega','_{\tau}',' (meV)'),'FontSize',16);
% xlabel('\omega_t (meV)','FontSize',16);
% axis equal;
% xlim([1641,1646]);
% ylim([-1646,-1641]);

% %build color map
% phaseamp = zeros(padi,padi,3);
% phaseamp(:,:,1) = ((atan2(imags1,reals1)+pi)/(2*pi)).*abs1/maxamp;
% phaseamp(:,:,3) = (1-atan2(imags1,reals1)-pi)/(2*pi).*abs1/maxamp;

figure();
%image(freq+refshift,freq-refshift,phaseamp);%,'edgecolor','none');%,15,'edgecolor','none');
imagesc(freq+refshift,freqtauref,atan2(imags1,reals1),'CDatamapping','scaled','AlphaData',(abs1/maxamp).^2);%1./(1+exp(-10*(abs1/maxamp-0.25))))
colormap('jet');
set(gca,'FontSize',14);
set(gca,'YDir','normal')
%title(strcat('Absolute of SI, P = ',Power,', T Delay = ',num2str(Tpos,2),', run ',num2str(Run)))
ylabel(strcat('\omega','_{\tau}',' (meV)'),'FontSize',16);
xlabel('\omega_t (meV)','FontSize',16);
axis equal;
xlim([freqref(1),freqref(padi)]);
% % % % xlim([1638.5,1643.5]);
% % % % ylim([-1643.5,-1638.5]);
% % % 
% % % % if ~s3
% % % % figure();
% % % % %image(freq+refshift,freq-refshift,phaseamp);%,'edgecolor','none');%,15,'edgecolor','none');
% % % % imagesc(freq+refshift,-freqtau+refshift,flipud(atan2(imags2,reals2)),'CDatamapping','scaled','AlphaData',(flipud(abs2/maxamp2)).^2);%1./(1+exp(-10*(abs2/maxamp2-0.45)))))
% % % % colormap('jet');
% % % % set(gca,'FontSize',14);
% % % % set(gca,'YDir','normal')
% % % % %title(strcat('Absolute of SI, P = ',Power,', T Delay = ',num2str(Tpos,2),', run ',num2str(Run)))
% % % % ylabel(strcat('\omega','_{\tau}',' (meV)'),'FontSize',16);
% % % % xlabel('\omega_t (meV)','FontSize',16);
% % % % axis equal;
% % % % colorbar;
% % % % xlim([freqref(1),freqref(padi)]);
% % % % end
% % % end
% % % 
% % % % % Plot X-diagonal slices
% % % % Ipk = [21,15,43]*4;
% % % % range = 10*4;
% % % % slices = zeros(2*range+1, 3);
% % % % sliceax = freq(1:2*range+1) - freq(range);
% % % % for i = 1:2*range+1
% % % %     slices(i,1) = reals1(i+(padi-Ipk(1))-range,i+Ipk(1)-range)/1.6730e4;
% % % %     slices(i,2) = reals1(i+(padi-Ipk(2))-range,i+Ipk(2)-range)/0.8227e4;
% % % %     slices(i,3) = reals1(i+(padi-Ipk(3))-range,i+Ipk(3)-range)/0.4989e4;
% % % % end
% % % % figure();
% % % % plot(sliceax,slices);
% % % % 
% % % % % Plot diagonal slice of SII
% % % % slices2 = zeros(padi,1);
% % % % for i = 1:padi
% % % %     slices2(i) = abs2(i,i);
% % % % end
% % % % figure();
% % % % plot(freq,slices2);
% % % % 
% % % % % figure();
% % % % % plot(freq,sum(abs1,2));
% % % 
% % % % Plot diagonal/X-diagonal
% % % [M,I] = max(abs1(:));
% % % [I_tau,I_t] = ind2sub(size(abs1),I);
% % % 
% % % homo = [1,1;padi,padi]; %[-tau,-t;+tau,+t]
% % % 
% % % maxi = homo(1,2)-I_t+I_tau; % for t limited
% % % if (maxi >= homo(1,1))
% % %     homo(1,1) = maxi;
% % % else % tau limited
% % %     homo(1,2) = homo(1,1)-I_tau+I_t;
% % % end
% % % maxi = homo(2,2)-I_t+I_tau;
% % % if (maxi <= homo(2,1))
% % %     homo(2,1) = maxi;
% % % else
% % %     homo(2,2) = homo(2,1)-I_tau+I_t;
% % % end
% % % 
% % % homoAx = zeros(abs(homo(2,1)-homo(1,1)+1), 1);
% % % for i = 1:(homo(2,1)-homo(1,1)+1)
% % %     homoAx(i) = abs1(i+homo(1,1)-1,i+homo(1,2)-1);
% % % end
% % % freqhomo = (freq(2)-freq(1))*linspace(-abs(homo(1,1)-I_tau),abs(homo(2,1)-I_tau),homo(2,1)-homo(1,1)+1);
% % % %(freq(2)-freq(1))*linspace(-abs((homo(1,1)-I_tau)+1j*(homo(1,2)-I_t)),abs((homo(2,1)-I_tau)+1j*(homo(2,2)-I_t)),homo(2,1)-homo(1,1)+1);
% % % 
% % % inhomo = [padi,1;1,padi];
% % % maxi = -(inhomo(1,2)-I_t)+I_tau; % for t limited
% % % if (maxi <= inhomo(1,1))
% % %     inhomo(1,1) = maxi;
% % % else % tau limited
% % %     inhomo(1,2) = -(inhomo(1,1)-I_tau)+I_t;
% % % end
% % % maxi = -(inhomo(2,2)-I_t)+I_tau;
% % % if (maxi >= inhomo(2,1))
% % %     inhomo(2,1) = maxi;
% % % else % tau limited
% % %     inhomo(2,2) = -(inhomo(2,1)-I_tau)+I_t;
% % % end
% % % 
% % % inhomoAx = zeros(abs(inhomo(2,1)-inhomo(1,1))+1, 1);
% % % for i = 1:abs(inhomo(2,1)-inhomo(1,1))+1
% % %     inhomoAx(i) = abs1(inhomo(1,1)-i+1,i+inhomo(1,2)-1);
% % % end
% % % freqinhomo = (freq(2)-freq(1))*linspace(-abs(inhomo(1,1)-I_tau),abs(inhomo(2,1)-I_tau),abs(inhomo(2,1)-inhomo(1,1))+1);
% % % %(freq(2)-freq(1))*linspace(-abs((inhomo(1,1)-I_tau)+1j*(inhomo(1,2)-I_t)),abs((inhomo(2,1)-I_tau)+1j*(inhomo(2,2)-I_t)),abs(inhomo(2,1)-inhomo(1,1))+1);
% % % figure();
% % % plot(freqhomo+freq(21)+refshift,homoAx);
% % % hold on;
% % % plot(freqinhomo+freq(21)+refshift,inhomoAx);
% % % title(strcat('S_I Slices, P = ',Power,', T Delay = ',num2str(Tpos,2),' ps, run ',num2str(Run)))
% % % set(gca,'FontSize',14);
% % % ylabel('S_I Amplitude (arb. units)')
% % % xlabel('\omega_t (meV)')
% % % %xlim([-5,5])
% % % 
% % % %fit function for (Lorentzian)
% % % lorFun = @(b,x)((b(1)*b(2))./((2*pi)*((x-b(3)).^2+(b(1)).^2)));
% % %     %b(1) = HWHM, b(2) = integral of amp, b(3) = peak pos
% % %     
% % % %Fit homo line
% % % beta0 = [0.1,1,0];
% % % beta0(2) = -sum(homoAx)*(freqhomo(2)-freqhomo(1));
% % % %W = squeeze(1./transpose(homoAx));
% % % 
% % % [betahomo,RT,JT,CovBT,MSET] = nlinfit(freqhomo,transpose(homoAx),lorFun,beta0);%,'Weights',W);
% % % CovBhomo = sqrt(diag(CovBT));
% % % hold on;
% % % plot(freqhomo+freq(21)+refshift,lorFun(betahomo,freqhomo));
% % % 
% % % %fit function for gaussian
% % % gaussasyFun = @(b,x)(b(2)*exp(-(x-b(3)).^2/(2*b(1)^2)).*(1+erf(b(4)*(x-b(3))/(sqrt(2)*b(1)))));
% % %     %b(1) = HWHM, b(2) = integral of amp, b(3) = peak pos, b(4) = skew
% % % lorasyFun = @(b,x)((b(1)*b(2))./((2*pi)*((x-b(3)).^2+(b(1)).^2)).*(1+erf(b(4)*(x-b(3))/(sqrt(2)*b(1)))));
% % %     
% % % %Fit inhomo line
% % % beta0 = [1,1,0,0];
% % % beta0(2) = max(inhomoAx);%-sum(inhomoAx)*(freqinhomo(2)-freqinhomo(1));
% % % %W = squeeze(1./transpose(inhomoAx));
% % % 
% % % [betainhomo,RT,JT,CovBT,MSET] = nlinfit(freqinhomo,transpose(inhomoAx),gaussasyFun,beta0);%,'Weights',W);
% % % CovBinhomo = sqrt(diag(CovBT));
% % % 
% % % %%%% Remove inhomo fit for now %%%%
% % % %hold on;
% % % %plot(freqinhomo,gaussasyFun(betainhomo,freqinhomo));
% % % 
% % % betahomo;
% % % betainhomo;
% % % homoW = betahomo(1);
% % % homoWe = CovBhomo(1);
% % % inhomoW = betainhomo(1)*sqrt(2*log(2));
% % % inhomoWe = CovBinhomo(1)*sqrt(2*log(2));
% % % wtau = freq(I_tau);
% % % wt = freq(I_t);
% % % save(strcat(folder,'Params',Run,'_',num2str(Tstep),'.mat'),'M','wtau','wt','homoW','homoWe','inhomoW','inhomoWe','Tpos');
% % % 
% % % % justt = abs(fftshift(squeeze(fft(procData(:,5,:)+1j*procData(:,7,:),padi,1))));
% % % % figure();
% % % % plot(freq,justt(:,1))