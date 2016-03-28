
function ClassDemo()
%yo()
main();
end

function main
clear all; close all; clc; dbstop if error;
sampleDir = {
    '2015-01-01_15-45-34_dreadd_rat_ref_animalground_2100depth'; %40 million samples
    '2015-01-29_dreadd_ref17_with_commutator_thresh50';          %15 million samples
    '2015-01-06_130_rat130_arena1' ;                             %50 million samples
    '2015-05-18_messi_before_injection_threshold40_ref4'};%dori  %60 million samples
file = char(sampleDir(4));
Fs = 32000;                   % Sampling frequency
T = 1/Fs;                     % Sample time
startCh = 1;                  % Start channel
numOfCh = 16;                 % End   channel
chunk = 30;                   % Chunk of 1e6 samples
numChunks = 2;                % Number of chunks

xlabFreq='Normalized frequency'; ylabFreq='Magnitude(\muV)';
xlabSig='millisecs'; ylabSig='\muV';

disp('start')


%loading from file (processed from raw data using NeuralDataBinaryToInteger.java script) 
loadFromFile=0;
%{
disp('ch loaded:');
for i = startCh:(startCh+numOfCh-1)
    channel(:,i) = loadChunk(file, i, chunk, numChunks);
    channel(:,i) = channel(:,i)*-1; %invert, makes spikes positive
    channel(:,i) = channel(:,i) -  mean(channel(:,i));
    fprintf('%d|',i);
end
%}
if loadFromFile == 0
    load('demo_channels');
    load('demo_fourier');
    load('demo_fourierCleanlNeg');
    load('demo_cleanChannelNeg');
    load('demo_HPfilteredCleanChannel');
    load('demo_spikes');
end
%}

% Time
time = (1:length(channel(:,1)))*T*1000; %converts to millisecs

%plot input
%all
figure; plotOffset(time,channel,500,'Input Signal: Raw Neural Recording of Rat',xlabSig, ylabSig);
%1 channel
figure; plotOffset(time,channel(:,9),500,'Channel 9: input Signal',xlabSig, ylabSig);


% Fourier
if loadFromFile == 1
    for i = startCh:(startCh+numOfCh-1)
        fourier(:,i) = fft(channel(:,i));
    
    end
end
disp('here');
%fourier cleaning (negatives)
for i = startCh:(startCh+numOfCh-1)
    tmp = cleanHarmonicsNeg(fourier(:,i));
    cleanFourierNeg(:,i) = tmp;
    cleanChannelNeg(:,i) = ifft(tmp);
end
tmp = 0;
disp('cleaned');

% plot fourier before and after after harmonics are cleaned
ftime = Fs/2*linspace(0,1,(length(fourier(:,1))/2+1));
tmp = abs(fourier(1:(length(fourier(:,1))/2)+1,1:length(fourier(1,:))));
set(gca,'DefaultTextFontSize',32);
figure; plot1xy(ftime, tmp(:,9), 'Fourier of input signal', xlabFreq, ylabFreq);
%plot both
tmp2 = abs(cleanFourierNeg(1:(length(cleanFourierNeg(:,1))/2)+1,1:length(cleanFourierNeg(1,:))));
figure; plotOverlayxy(ftime,tmp(:,9),tmp2(:,9), 'Fourier after harmonics cleaned', xlabFreq, ylabFreq);
legend('original', 'after cleaning harmonics');

% plot channel after its been cleaned of harmonics
figure; plotOffset(time,cleanChannelNeg,500,'Harmonics cleaned channels',xlabSig, ylabSig);
%just channel 9
figure; plotOffset(time,cleanChannelNeg(:,9),500,'Channel 9: harmonics cleaned',xlabSig, ylabSig);
%original vs cleanes
figure; plotOverlayxy(time,channel(:,9), real(cleanChannelNeg(:,9)),'Channel 9: harmonics cleaned',xlabSig, ylabSig);
legend('original', 'after cleaning harmonics');


%%% High Pass 
%hp only, for debug
if loadFromFile == 1
    for numSpikesPerCluster = startCh:(startCh+numOfCh-1)
        %HPfilteredChannel(:,i) = highPassFilter(channel(:,i), Fs);
    end
end

% hp filter plus cleaned
if loadFromFile == 1
    for i = startCh:(startCh+numOfCh-1)
        HPfilteredCleanChannel(:,i) = highPassFilter(cleanChannelNeg(:,i), Fs);
    end
end
disp('filtered');
%plotting fourier original, cleaning harmonics, and highpass filtering.
tmp3=fft(HPfilteredCleanChannel(:,9));
tmp4 = abs(tmp3(1:(length(tmp3)/2)+1));
figure;
hold on
plot1xy(ftime,tmp(:,9), '', xlabFreq, ylabFreq);
plot1xy(ftime,tmp2(:,9), '', xlabFreq, ylabFreq);
plot1xy(ftime,tmp4, 'Fourier after Highpass Filtering', xlabFreq, ylabFreq);
hold off
legend('original', 'after cleaning harmonics', 'after highpass filter of 600Hz');
tmp = 0; tmp2 = 0, tmp3 = 0, tmp4 = 0;

% plot High Pass
figure; plotOffset(time,HPfilteredCleanChannel,500,'Harmonics cleaned and high-pass filtered channels',xlabSig, ylabSig);
% Channel 9
figure; plot1xy(time,HPfilteredCleanChannel(:,9),'Channel 9: harmonics cleaned and filtered', xlabSig, ylabSig);
% overlay
figure; plotOverlayxy(time,real(channel(:,9)), real(HPfilteredCleanChannel(:,9)),'Channel 9: after harmonics cleaned and filtered',xlabSig, ylabSig);
legend('original', 'after cleaning and highpass filter');

%%%%%%%%%       Extracting spikes 
if loadFromFile == 1
    spikes = extractSpikes(HPfilteredCleanChannel, 70, 32, 64);
end
disp('spiked extracted');
nsp = zeros(1,16);
for numSpikesPerCluster = 1:length(spikes(2,:))
    nsp(spikes(2,numSpikesPerCluster)) = nsp(spikes(2,numSpikesPerCluster)) + 1;
end
%plot spikes per channel
figure; bar(nsp./sum(nsp)*100), title('% spikes per channel'); % which channel had the most spikes
%cleaning high values
vthreshold = 400;
signals = spikes(3:end,:);
signals(signals>vthreshold)=[vthreshold];
signals(signals<-vthreshold)=[-vthreshold];

%kmeans input
nki = 6;
tic
kmi = kmeans(signals',nki);
disp('elapsed time K-means on input');
tIn = toc
[n , p] = hist(kmi,nki); bar(1:nki,n./sum(n)*100), title(sprintf('kmeans on events k=%d non-svd', nki));
ylabel('% of total events'); xlabel('cluster id');

%%% sort signals by cluster
clusters = cell(nki,16);
numSpikesPerCluster=zeros(nki,1);
sampleLength = 96;
for k = 1:length(kmi)
  numSpikesPerCluster(kmi(k)) = numSpikesPerCluster(kmi(k))+1;
  for ch = 1:16 %16 channels
  	clusters{kmi(k),ch}(end+1,:) = signals(1+(ch-1)*sampleLength:ch*sampleLength,k);
  end
end
%%%

%%%plot input
%%%mean
close all;
for k = 1:size(clusters,1)
    figure(k);
    for ch = 1:size(clusters,2)
        if length(clusters{k,1}) > 100
        subplot(4,4,ch)
        for sam = 1:sampleLength
            cmean(k,ch,1,sam) = mean(clusters{k,ch}(:,sam));
            cmean(k,ch,2,sam) =  std(clusters{k,ch}(:,sam));
        end
        errorbar(cmean(k,ch,1,:), cmean(k,ch,2,:)),
        axis([1,98,-150,150]);
        if k==1
            %title(sprintf('k_id:%d ch:%d n:%d non-svd', k, ch,numSpikesPerCluster(k)));
        else
            %title('');
        end
        
        end
    end
    suptitle(sprintf(' non-svd: kluster:%d channel:%d n-spikes:%d\n y=microvolts x=time',...
        k, ch,numSpikesPerCluster(k)));
    saveas(gcf,sprintf('demo_ns_mean_k_%d.png',k));
end

%Single Value Decomposition %SVD%
%make sure data matrix is set up so samples are COL
% u is basis vector
%then multiply first n COL of u TRANSPOSE times datam atrix to change base 
%do svd on that.
[u , s , v] = svd(signals,'econ');
figure; plot(s), title('Singular values of data SVD');
nUvectors = 40;
pcaSpikes = u(:,1:nUvectors)'*signals;%converts spikes to new base using only u coefficients
%original vectors
originalBasisSpikes = u(:,1:nUvectors)* pcaSpikes;
mse = mean(mean((signals(:,:)- originalBasisSpikes(:,:)).^2));
%memory storage
figure; bar([100,numel(pcaSpikes)/numel(signals)*100]), 
title(sprintf('Storage of original signal vs. SVD using %d u vectors',nUvectors));
set(gca,'XTick',1:2,'XTickLabel',{'original';'pca'}), ylabel('% of original');
%figure original vs pca
figure;
hold on;
plot(signals(:,700));
plot(originalBasisSpikes(:,700));
hold off;
title(sprintf('Sample of original vs SVD Vector using %d u vectors\n(average MSE for all samples %.2f)',nUvectors,mse));
legend('original', 'svd');
%sanity check MSE per # of u's 
%mseplot = zeros(2,length(u)/10+1);
for ui=1:10:length(u)
    pcaSpikes2 = u(:,1:ui)'*signals;
    originalBasisSpikes2 = u(:,1:ui)* pcaSpikes2;
    uii = ui;
    if(ui > 1)
        uii = (ui-1)/10 +1;
    end
    disp(ui);
    mseplot(1,uii) = mean(mean((signals(:,:)- originalBasisSpikes2(:,:)).^2));
    mseplot(2,uii) = ui;
end
figure;plot(mseplot(2,:),mseplot(1,:)); title('MSE per number of SVD u vectors');
ylabel('mse'), xlabel('number of u vectors');

%%%K%%%%%%%%%%%%%%%%%% K-Means after PCA
nk = 6;
tic
kmPca = kmeans(pcaSpikes',nk);
disp('elapsed time K-means on pca');
tPca = toc;
%time:
%Elapsed time is 8.505078 seconds.
%Elapsed time is 0.110036 seconds.
figure; bar([tIn,tPca]), 
title(sprintf('Time of K-Means on original vs SVD using %d u vectors',nUvectors));
set(gca,'XTick',1:2,'XTickLabel',{'original';'pca'}), ylabel('sec');
%spikes per cluster pca
figure; hist(kmPca), title('number of kmeans clusters pca signals');
%n = [342, 308, 348, 292, 10, 303]; nk = 6;
figure; [n , p] = hist(kmPca,nk); bar((1:nk),n./sum(n)*100), ...
    title(sprintf('kmeans on events k=%d SVD', nk));
ylabel('% of total events'); xlabel('cluster id');


%%% sort PCA signals by cluster
clustersPCA = cell(nk,16);
numSpikesPerClusterPCA=zeros(nk,1);
sampleLength = 96;
for k = 1:length(kmPca)
  numSpikesPerClusterPCA(kmPca(k)) = numSpikesPerClusterPCA(kmPca(k))+1;
  for ch = 1:16 %16 channels
  	clustersPCA{kmPca(k),ch}(end+1,:) = signals(1+(ch-1)*sampleLength:ch*sampleLength,k);
  end
end
%%%

%%%plot PCA ckusters
%%%mean
close all;
for k = 1:size(clustersPCA,1)
    figure(k);
    for ch = 1:size(clustersPCA,2)
        if length(clustersPCA{k,1}) > 100
        subplot(4,4,ch)
        for sam = 1:sampleLength
            cmean(k,ch,1,sam) = mean(clustersPCA{k,ch}(:,sam));
            cmean(k,ch,2,sam) =  std(clustersPCA{k,ch}(:,sam));
        end
        errorbar(cmean(k,ch,1,:), cmean(k,ch,2,:)),
        axis([1,98,-150,150]);
        if k==1
            %title(sprintf('k_id:%d ch:%d n:%d non-svd', k, ch,numSpikesPerCluster(k)));
        else
            %title('');
        end
        
        end
    end
    suptitle(sprintf(' SVD: kluster:%d channel:%d n-spikes:%d\n y=microvolts x=time',...
        k, ch,numSpikesPerClusterPCA(k)));
    saveas(gcf,sprintf('demo_SVD_mean_k_%d.png',k));
end


stop
stop
end


%%%%%%%%%%%SUBROUTINES%%%%%%%%%%
function sortSpikesByClusterChan(ki, signals)
    sampleLength = 96;
    for k = 1:length(ki)
        for ch = 1:16 %16 channels 
            sorted(ki(k),ch,1:sampleLength) = signals(1+(ch-1)*sampleLength:ch*sampleLength,k); 
        end
    end
end


%MAKE SURE DOESNT FAST FORWARD OTHER CHANNELS
function spikes = extractSpikes(data, thresholdMiV, windowBeforeMS, windowAfterMS)
%data: each channel is column vectors
%start at first ms so window doesn't crash
%data matrix, rows = sample structure, cols = samples
%spikes = double(zeros(2+(windowBeforeMS+windowAfterMS)*16,countSpikes(data,thresholdMiV, windowBeforeMS, windowAfterMS)));
eventIntervalThresh = 32;%1ms between events
spikes = zeros(2+(windowBeforeMS+windowAfterMS)*16,10e3); %max number make larger
numSpikes = 2; %lameness
for timei = windowBeforeMS+1:length(data(:,1))-windowAfterMS
   saved = 0;
   %channel loop
   for channeli = 1: length(data(1,:)) %could optimize %&& (saved < 1)
       if (data(timei,channeli) > thresholdMiV) && ... 
       (timei > (spikes(1,numSpikes-1)+ eventIntervalThresh)) %space between event intervals
          saved = 1;
          %fprintf('%d > %d :%d\n',timei, spikes(1,numSpikes)+ windowAfterMS,(timei > (spikes(1,numSpikes)+ windowAfterMS)));
          spikes(1:2,numSpikes) = [timei;channeli]; %store time and channel of spike recorded
          for chi2 = 1: length(data(1,:))
              %crazy indexing! +1 is to get next index +2 is the offset of
              %storing time and channel is first two values
              spikes((1+2+(chi2-1)*(windowBeforeMS+windowAfterMS)):2+chi2*(windowBeforeMS+windowAfterMS),numSpikes) = ...
                  real(data((timei-windowBeforeMS):(timei+windowAfterMS-1),chi2)); %32 samples per ms
          end
          numSpikes = numSpikes + 1;
       %else
          % disp('skipped');
       end
   end
end
spikes = spikes(:,2:numSpikes-1);
end

% Clean Harmonics
%this function cleans harmonics by picking a starting point(offset) after
%which any value above a threashold (mean*factor) is reset to mean.
function out=cleanHarmonicsNeg(input)
offset = 2*1e4; %when to start filtering, this is determined experimentally 
% looking at the original fourier and seeing that harmonics do not start until these frequencies
tmp=abs(input); %to determine the harmonics taking the abs is necessary
out=input;
factor = 3; %factor times mean to set threshold for cutting (determined experimentally)
period = 100;%the period over which to determine the average 
sumPer = sum(tmp(offset+1:offset+period)); 
for i = offset+period+1:length(tmp) 
    if tmp(i)>factor*sumPer/period %mean over period
        tmp(i)=sumPer/period;  %changes value to mean
        out(i)= sumPer/period;
        if out(i)<0
            out(i)= -1*sumPer/period;
        end
    end
    sumPer=sumPer-tmp(i-period)+tmp(i);
end
end

function data = loadChannel(file, channel)
hdir = 'C:\\Users\\alm\\Desktop\\dori\\raw_data';
ADBitVolts = 0.000000015624999960550667;
data = load(sprintf('%s\\%s\\ch%d.csv',hdir,file,channel));
data = data*ADBitVolts; % to volts
data = data*1e6; % to micro volts
end

%chunks are 1e6 long consequetive parts of the channel signal, can load
%multiple chunks
function data = loadChunk(file, channel, chunk, numChunks)
hdir = 'C:\\Users\\alm\\Desktop\\dori\\raw_data';
ADBitVolts = 0.000000015624999960550667;
data = zeros(1e6*numChunks,1); %chunks are 1e6 long, preallocates
for i = 1:numChunks
    tmp = load(sprintf('%s\\%s\\chunks\\ch%d_%d.csv',hdir,file,channel,chunk+i-1));
    data((i-1)*1e6+1:i*1e6) = tmp(:,2);
end
data = data*ADBitVolts; % to volts % to volts 
data = data*1e6; % to micro volts
%data = data';
end

function data = highPassFilter(input, samplingFrequency)
%HIGH PASS FILTERING
% http://www.mathworks.com/help/dsp/ref/fdesign.bandpass.html
% All frequency values are in Hz.
% Construct an FDESIGN object and call its BUTTER method.
Fstop =590;           % First Stopband Frequency
Fpass =600;           % First Passband Frequency
Astop = 10;          % First Stopband Attenuation (dB)
Apass  = 5;           % Passband Ripple (dB)
fs = samplingFrequency;
%toDecibal = 20*log(10);         
Hd = design(fdesign.highpass(Fstop, Fpass, Astop, Apass, fs),'butter');
data = filtfilt(Hd.sosMatrix,Hd.ScaleValues,input);
end

%%%%%%%%%PLOTTING FUNCTIONS%%%%%%%%%%



function plot1xy(x,y,tit, xlab, ylab)
    plot(x,y)
    title(tit)
    xlabel(xlab)
    ylabel(ylab)
end

function plot1(y, tit, xlab, ylab)
    plot(y)
    title(tit)
    xlabel(xlab)
    ylabel(ylab)
end

%for plotting frequency
function plotOverlay(data, data2, tit, xlab, ylab, cols)

nplots = size(data,2);
if nplots > 1
    for i = 1:nplots %assumes 4 plots/adata
        subplot(nplots/cols,cols,i);
        hold on
        plot(data(:,i))
        plot(data2(:,i))
        hold off
        title(tit)
        xlabel(xlab)
        ylabel(ylab)
    end
else
     hold on
     plot(data)
     plot(data2)
     hold off
end
title(tit)
xlabel(xlab)
ylabel(ylab)
end

%for plotting frequency
function plotOverlayxy(x,data, data2, tit, xlab, ylab, cols)

nplots = size(data,2);
if nplots > 1
    for i = 1:nplots 
        subplot(nplots/cols,cols,i);
        hold on
        plot(x,data(:,i))
        plot(x,data2(:,i))
        hold off
        title(tit)
        xlabel(xlab)
        ylabel(ylab)
    end
else
     hold on
     plot(x,data);
     plot(x,data2);
     hold off
end
title(tit)
xlabel(xlab)
ylabel(ylab)
end



%for plotting frequency
function plotSubplots(data, tit, xlab, ylab, cols)
nplots = size(data,2);
    for i = 1:nplots %assumes 4 plots/adata
        subplot(nplots/cols,cols,i);
        plot(data(:,i))
        title(tit)
        xlabel(xlab)
        ylabel(ylab)
    end
end

%for plotting signal, plots all sets vertically offset
function plotOffset(time, data, offset, tit, xlab, ylab)
    hold on;
    for i = length(data(1,:)):-1:1
         plot(time,data(:,i)+(i-1)*offset)
         %plot(time,ones(1,length(data(:,i)))*(i-1)*offset)
    end
    hold off
    title(tit);
    xlabel(xlab);
    ylabel(ylab);
    legend('show');
    %n=get(gca,'Ytick');
    %set(gca,'yticklabel',sprintf('%.0f',n'));
end


%for plotting frequency
function plot4subplotsOverlay(data, data2, tit, xlab, ylab)
    for i = 1:4 %assumes 4 plots/adata
        subplot(2,2,i);
        hold on
        plot(data(:,i))
        plot(data2(:,i))
        hold off
        title(tit)
        xlabel(xlab)
        ylabel(ylab)
    end
end

%for plotting frequency
function plot4subplots(data, tit, xlab, ylab)
    for i = 1:4 %assumes 4 plots/adata
        subplot(2,2,i);
        plot(data(:,i))
        title(tit)
        xlabel(xlab)
        ylabel(ylab)
    end
end
