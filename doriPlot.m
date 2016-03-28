
function doriPlot()
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
chunk = 1;                    % Chunk of 1e6 samples
numChunks = 3;                % Number of chunks

xlabFreq='Normalized frequency'; ylabFreq='Magnitude(\muV)';
xlabSig='millisecs'; ylabSig='\muV';

disp('start')
disp('ch loaded:');
%loading
for i = startCh:(startCh+numOfCh-1)
    channel(:,i) = loadChunk(file, i, chunk, numChunks);
    channel(:,i) = channel(:,i)*-1; %invert, makes spikes positive
    fprintf('%d|',i);
end

% Time
time = (1:length(channel(:,1)))*T*1000; %converts to millisecs

% dori vectors
%channelPlot(:,startCh:startCh+numOfCh-1) = ...
%    channel(:,startCh:startCh+numOfCh-1) + offset;
%figure; plot(time,channel); title('input signal');

%plot input
%figure; plotOffset(time,channel,300,'input signal',xlabSig, ylabSig);

% Fourier
for i = startCh:(startCh+numOfCh-1)
    tmp = fft(channel(:,i));
    fourier(:,i) = tmp;
    %data(:,1) = (abs(tmp(1:(length(tmp)/2)+1)));
end
% plot fourier of original signal 
%figure; plot4subplots(data, 'absolute value fourier of input', xlabFreq, ylabFreq);

%fourier cleaning (negatives)
for i = startCh:(startCh+numOfCh-1)
    %don't need %tmp2 = cleanHarmonics(fourier(:,i));%cleanFourierAbs(:,i) = tmp2;%cleanChannelAbs(:,i) = ifft(tmp2);
    tmp = cleanHarmonicsNeg(fourier(:,i));
    %cleanFourierNeg(:,i) = tmp;
    cleanChannelNeg(:,i) = ifft(tmp);
    %data(:,i) = (abs(tmp(1:(length(tmp)/2)+1)));     
end
disp('cleaned');
fourier = 0; tmp = 0;
% plot fourier before and after after harmonics are cleaned
%tmp = abs(fourier(1:(length(fourier(:,1))/2)+1,1:length(fourier(1,:))));
%figure; plot4subplots(tmp, 'Fourier: cleaned harmonics', xlabFreq, ylabFreq);
       %figure; plot4subplots(data,'Fourier: unclean vs cleaned harmonics', xlabFreq, ylabFreq);

% plot channel after its been cleaned
%figure; plot4spaced(time,cleanChannelNeg,'cleaned inversed fourier channels',xlabSig, ylabSig)

%%% High Pass filter only
for i = startCh:(startCh+numOfCh-1)
    %HPfilteredChannel(:,i) = highPassFilter(channel(:,i), Fs);
end
% hp filter plus cleaned
for i = startCh:(startCh+numOfCh-1)
    HPfilteredCleanChannel(:,i) = highPassFilter(cleanChannelNeg(:,i), Fs);
end
disp('filtered');
cleanChannelNeg = 0;

%%%  Band Pass filter
for i = startCh:(startCh+numOfCh-1)
    %bandFilteredCleanChannel(:,i) = bandPassFilter(cleanChannelNeg(:,i), Fs);
end
% plot High Pass
%figure; plot4spaced(time,HPfilteredCleanChannel,'clean and high-pass filtered channels',xlabSig, ylabSig);
% plot Band filtered
%figure; plot4spaced(time,bandFilteredCleanChannel,'Band filtered channels',xlabSig, ylabSig);
% both
%figure; plot4spaced(time,HPfilteredCleanChannel,[],[],[]);
%plot4spaced(time,bandFilteredCleanChannel,'High and Band filtered channels',xlabSig, ylabSig);

%%% Weiner
for i = startCh:(startCh+numOfCh-1)
    %sig = HPfilteredCleanChannel(:,i); 
    %sigma = 0.03*(max(sig)-min(sig));
    %wiener(:,i) = WienerFilter(sig,sig,sigma);
end
% plot Weiner and High Pass
%figure; plot4spaced(time,wiener,'wiener filtered channels',xlabSig, ylabSig);
% before after wiener
%figure; plot4spaced(time,HPfilteredCleanChannel,[],[],[]);
%plot4spaced(time,wiener,'before/after wiener filter (from highpass)',xlabSig, ylabSig);

tmp = 0; fourier = tmp; cleanChannelNeg = tmp; channel = tmp;
%figure; plotOffset(time,HPfilteredCleanChannel, 300, sprintf('%s:\nclean and high-pass filtered channels', file),xlabSig, ylabSig); 
%add name to variables
spikes = extractSpikes(HPfilteredCleanChannel, 70, 32, 64);
disp('spiked extracted');
[n, p] = hist(spikes(2,:),16);
figure; bar(p,n./sum(n)*100), title('spikes per channel'); % which channel had the most spikes

%Single Value Decomposition %SVD%
HPfilteredCleanChannel = [];
vthreshold = 400;
signals = spikes(3:end,:);
signals(signals>vthreshold)=[vthreshold];
signals(signals<-vthreshold)=[-vthreshold];

%make sure data matrix is set up so samples are COL
% u is basis vector
%then multiply first n COL of u TRANSPOSE and do svd on that.
[u , s , v] = svd(signals,'econ');
figure; plot(s), title('singular values of signals');
%principleComponents = u(:,1:50)*s(1:50,1:50);
nk = 10;
pcaSpikes = u(:,1:100)'*signals;%converts spikes to new base
kmPca = kmeans(pcaSpikes',nk);
figure; hist(kmPca), title('kmeans pca signals');
originalBasisSpikes = u(:,1:100)* pcaSpikes;

figure;
hold on;
plot(signals(1000,:));
%figure
plot(originalBasisSpikes(1000,:));
hold off;


%u, k loop
for ui = 150%:50:300
   figure;
   c = 1;
   for ki = 10:1:14
       pcaSpikes = u(:,1:ui)'*signals;
       kmPca = kmeans(pcaSpikes',ki);
       subplot(5,1,c); c = c+1;
       [n , p] = hist(kmPca);
       bar(p,n./sum(n)*100), title(sprintf('u:%d k:%d', ui, ki));
   end
end

%kmeans input SLOW
nki = 10;
kmi = kmeans(signals',nki);
[n , p] = hist(kmi,nki);
bar(p,n./sum(n)*100), title(sprintf('kmeans on events k=%d', nki));
ylabel('% of total events');

%figure;hist(kmi), title('kmeans input');
%%% sort signals by cluster
clusters = cell(nki,16);
i=zeros(nki,1);
sampleLength = 96;
for k = 1:9546
  i(kmi(k)) = i(kmi(k))+1;
  for ch = 1:16 %16 channels
  	clusters{kmi(k),ch}(end+1,:) = signals(1+(ch-1)*sampleLength:ch*sampleLength,k);
  end
end
%%%

%spikes moving through time is an issue

%%%plot
close all;
for k = 1:size(clusters,1)
    figure(k);
    for ch = 1:size(clusters,2)
        fprintf('k:%d ch:%d\n', k, ch);
        subplot(4,4,ch)
        hold on
        for sam = 1:1:size(clusters{k,ch},1)
            %disp(sam); 
            %plot(clusters{k,ch,sam,:}), title(sprintf('k:%d ch:%d', k, ch));
            axis([1,98,-200,200]);
            plot(clusters{k,ch}(sam,:)), title(sprintf('kl:%d ch:%d sz:%d', k, ch,i(k)));
            saveas(gcf,sprintf('all_signals_k_%d.png',k));
        end
        hold off
    end
end

%%%mean
close all;
for k = 1:size(clusters,1)
    figure(k);
    for ch = 1:size(clusters,2)
        subplot(4,4,ch)
        for sam = 1:sampleLength
            cmean(k,ch,1,sam) = mean(clusters{k,ch}(:,sam));
            cmean(k,ch,2,sam) =  std(clusters{k,ch}(:,sam));
        end
        errorbar(cmean(k,ch,1,:), cmean(k,ch,2,:)),
        axis([1,98,-150,150]);
        title(sprintf('kl:%d ch:%d sz:%d', k, ch,i(k)));
        saveas(gcf,sprintf('error_bar_k_%d.png',k));
    end
end



hold on
% plots signals overlayed 
% figure;
for i = 1:10:9546
    %plot(signals(i,:));
end
hold off;

stop
stop
end

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

function spks = countSpikes(data, thresholdMiV, windowBeforeMS, windowAfterMS)
   spks = 0;
   for i = 1:length(data(1,:))
       spks = spks + sum(data(windowBeforeMS+1:length(data(:,1))-windowAfterMS,i)>thresholdMiV);
   end
   samples = 16*(3*32); %32 samples per ms 1 ms before, 2 after * 16 channels
   fprintf('%d spikes X %d samples %d elemets\n',spks,samples,samples*spks);
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

function data = bandPassFilter(input, samplingFrequency)
% http://www.mathworks.com/help/dsp/ref/fdesign.bandpass.html
% All frequency values are in Hz.
% Construct an FDESIGN object and call its BUTTER method.
Fstop =590;           % First Stopband Frequency
Fpass =600;           % First Passband Frequency
Astop = 20*log10(200);          %10 First Stopband Attenuation (dB)
Apass  = 20*log10(100);             %5 Passband Ripple (dB)

F2pass =5e3-50; %BE CAREFUL OF HIGH FREQUENCIES
F2stop =5e3;           % First Stopband Frequency

fs = samplingFrequency;
%toDecibal = 20*log(10);         
%d =        fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,C',100,800,1e3,1.4e3,1.6e3,1e4);
Bd = design(fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', Fstop, Fpass,F2pass, F2stop, Astop, Apass, Astop, fs),'butter');
data = filtfilt(Bd.sosMatrix,Bd.ScaleValues,input);
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

%for plotting signal, plots all sets vertically offset
function plotOffset(time, data, offset, tit, xlab, ylab)
    hold on;
    for i = length(data(1,:)):-1:1
         plot(time,data(:,i)+(i-1)*offset)
         plot(time,ones(1,length(data(:,i)))*(i-1)*offset)
    end
    hold off
    title(tit);
    xlabel(xlab);
    ylabel(ylab);
    legend('show');
    %n=get(gca,'Ytick');
    %set(gca,'yticklabel',sprintf('%.0f',n'));
end

%%%%%%%%%%%%%%%%NOT USED%%%%%%%%%%%%%%%%%

function yo()
bandPassFilter((1:1e7), 32000);
disp('done yo');
end


function out=cleanHarmonics(input)
offset = 2*1e4;
out=abs(input);%input(1:length(input)/2+1);%(offset:end);  %REMOVE!
period = 100;
%mean100 = mean(out(1:100));
sumPer = sum(out(offset+1:offset+period)); 
for i = offset+period+1:length(out) %2005050
    if out(i)>3*sumPer/period       %20000
        out(i)=sumPer/period;
    end
    sumPer=sumPer-out(i-period)+out(i);
end
end


function spikes = extractSpikesOrig(data, thresholdMiV, windowBeforeMS, windowAfterMS)
%data: each channel is column vectors
%start at first ms so window doesn't crash
%data matrix, rows = sample structure, cols = samples
spikes = double(zeros(2+(windowBeforeMS+windowAfterMS)*16,countSpikes(data,thresholdMiV, windowBeforeMS, windowAfterMS)));
numSpikes = 1;
for timei = windowBeforeMS+1:length(data(:,1))-windowAfterMS
   saved = 0;
   %channel loop
   for channeli = 1: length(data(1,:)) %could optimize
       if (data(timei,channeli) > thresholdMiV) %&& (saved < 1)
          saved = 1;
          spikes(1:2,numSpikes) = [timei;channeli]; %store time and channel of spike recorded
          for chi2 = 1: length(data(1,:))
              %crazy indexing! +1 is to get next index +2 is the offset of
              %storing time and channel is first two values
              spikes((1+2+(chi2-1)*(windowBeforeMS+windowAfterMS)):2+chi2*(windowBeforeMS+windowAfterMS),numSpikes) = ...
                  real(data((timei-windowBeforeMS):(timei+windowAfterMS-1),chi2)); %32 samples per ms
          end
          numSpikes = numSpikes + 1;
       end
   end
end

end

function spikes = extractSpikesNoTrace(data, thresholdMiV, windowBeforeMS, windowAfterMS)
spikes = double(zeros(2+(windowBeforeMS+windowAfterMS)*16,countSpikes(data,thresholdMiV, windowBeforeMS, windowAfterMS)));
numSpikes = 1;
timei = windowBeforeMS+1
%for timei = windowBeforeMS+1:length(data(:,1))-windowAfterMS
while timei < length(data(:,1))-windowAfterMS
   saved = 0;
   %channel loop
   for channeli = 1: length(data(1,:)) %could optimize
       if (data(timei,channeli) > thresholdMiV) && (saved < 1)
          saved = 1;
          spikes(1:2,numSpikes) = [timei;channeli]; %store time and channel of spike recorded
          for chi2 = 1: length(data(1,:))
              %crazy indexing! +1 is to get next index +2 is the offset of
              %storing time and channel is first two values
              spikes((1+2+(chi2-1)*(windowBeforeMS+windowAfterMS)):2+chi2*(windowBeforeMS+windowAfterMS),numSpikes) = ...
                  real(data((timei-windowBeforeMS):(timei+windowAfterMS-1),chi2)); %32 samples per ms
          end
          numSpikes = numSpikes + 1;   
       end
   end
   if saved > 0 %gets rid of traces
       timei = timei + windowAfterMS; %only save one spike per window
   else
       timei = timei +1;
   end
end
end
