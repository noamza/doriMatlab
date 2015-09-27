
function doriPlot()
%yo()
main();
end

function yo()
bandPassFilter((1:1e7), 32000);
disp('done yo');
end

function main
clear all; close all; clc;
sampleDir = {
    '2015-01-01_15-45-34_dreadd_rat_ref_animalground_2100depth'; %40 million samples
    '2015-01-29_dreadd_ref17_with_commutator_thresh50';          %15 million samples
    '2015-01-06_130_rat130_arena1' ;                             %50 million samples
    '2015-05-18_messi_before_injection_threshold40_ref4'};       %60 million samples
file = char(sampleDir(1));
Fs = 32000;                   % Sampling frequency
T = 1/Fs;                     % Sample time
startCh = 1;
numOfCh = 4;

xlabFreq='Normalized frequency'; ylabFreq='Magnitude(\muV)';
xlabSig='millisecs'; ylabSig='\muV';

%loading
for i = startCh:(startCh+numOfCh-1)
    channel(:,i) = loadChunk(file, i, 1);
end

% Time
time = (1:length(channel(:,1)))*T*1000;

% dori vectors
%channelPlot(:,startCh:startCh+numOfCh-1) = ...
%    channel(:,startCh:startCh+numOfCh-1) + offset;
%figure; plot(time,channel); title('input signal');

figure; plot4spaced(time,channel,'input signal',xlabSig, ylabSig);

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
    cleanFourierNeg(:,i) = tmp;
    cleanChannelNeg(:,i) = ifft(tmp);
    data(:,i) = (abs(tmp(1:(length(tmp)/2)+1)));     
end
% plot fourier before and after after harmonics are cleaned
tmp = abs(fourier(1:(length(fourier(:,1))/2)+1,1:length(fourier(1,:))));
hold on;
figure; plot4subplots(tmp, 'Fourier: cleaned harmonics', xlabFreq, ylabFreq);
        plot4subplots(data,'Fourier: unclean vs cleaned harmonics', xlabFreq, ylabFreq);
hold off;

% plot channel after its been cleaned
figure; plot4spaced(time,cleanChannelNeg,'cleaned inversed fourier channels',xlabSig, ylabSig)

%%% High Pass filter only
for i = startCh:(startCh+numOfCh-1)
    HPfilteredChannel(:,i) = highPassFilter(channel(:,i), Fs);
end
% hp filter plus cleaned
for i = startCh:(startCh+numOfCh-1)
    HPfilteredCleanChannel(:,i) = highPassFilter(cleanChannelNeg(:,i), Fs);
end

%%%  Band Pass filter
for i = startCh:(startCh+numOfCh-1)
    bandFilteredCleanChannel(:,i) = bandPassFilter(cleanChannelNeg(:,i), Fs);
end
% plot High Pass
figure; plot4spaced(time,HPfilteredCleanChannel,'clean and high-pass filtered channels',xlabSig, ylabSig);
% plot Band filtered
figure; plot4spaced(time,bandFilteredCleanChannel,'Band filtered channels',xlabSig, ylabSig);
% both
figure; plot4spaced(time,HPfilteredCleanChannel,[],[],[]);
plot4spaced(time,bandFilteredCleanChannel,'High and Band filtered channels',xlabSig, ylabSig);

%%% Weiner
for i = startCh:(startCh+numOfCh-1)
    sig = HPfilteredCleanChannel(:,i); 
    sigma = 0.03*(max(sig)-min(sig));
    wiener(:,i) = WienerFilter(sig,sig,sigma);
end
% plot Weiner and High Pass
figure; plot4spaced(time,wiener,'wiener filtered channels',xlabSig, ylabSig);
% before after wiener
figure; plot4spaced(time,HPfilteredCleanChannel,[],[],[]);
plot4spaced(time,wiener,'before/after wiener filter (from highpass)',xlabSig, ylabSig);

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

%chunks are 1e6 long consequetive parts of the channel signal
function data = loadChunk(file, channel, chunk)
hdir = 'C:\\Users\\alm\\Desktop\\dori\\raw_data';
ADBitVolts = 0.000000015624999960550667;
data = load(sprintf('%s\\%s\\chunks\\ch%d_%d.csv',hdir,file,channel,chunk));
data = data(:,2)*(ADBitVolts*1e6); % to volts % to micro volts
%data(:,1) = data(:,2)*1e6; 
data = data';
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

%for plotting signal
function plot4spaced(time, data, tit, xlab, ylab)
    offset = 5000; %y axis MAGIC NUMBER
    hold on;
    for i = 1:4
         plot(time,data(:,i)+i*5000)
    end
    hold off
    title(tit);
    xlabel(xlab);
    ylabel(ylab);
    legend('show');
    %n=get(gca,'Ytick');
    %set(gca,'yticklabel',sprintf('%.0f',n'));
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

