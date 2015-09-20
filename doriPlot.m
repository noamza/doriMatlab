
function doriPlot()
%yo()
main();
end

function yo()
cleanHarmonics((1:1e7));
end

function main
clear all; close all; clc;
%{
2015-01-01_15-45-34_dreadd_rat_ref_animalground_2100depth     40 million samples
2015-01-29_dreadd_ref17_with_commutator_thresh50              15 million samples
2015-01-06_130_rat130_arena1                                  50 million samples
2015-05-18_messi_before_injection_threshold40_ref4            60 million samples
%}
sampleDir = {'2015-01-01_15-45-34_dreadd_rat_ref_animalground_2100depth'; 
         '2015-01-29_dreadd_ref17_with_commutator_thresh50'; 
         '2015-01-06_130_rat130_arena1' ;
         '2015-05-18_messi_before_injection_threshold40_ref4'};
file = char(sampleDir(1));
Fs = 32000;                   % Sampling frequency
T = 1/Fs;                     % Sample time
startCh = 1;
numOfCh = 4;

channel = zeros(1e6,4);
t = loadChunk(file, 1, 1);
channel(:,1)=t;
%loading
for i = startCh:(startCh+numOfCh-1)
    channel(:,i) = loadChunk(file, i, 1);
end

offset = 5000;

%channelPlot(:,startCh:startCh+numOfCh-1) = ...
 %   channel(:,startCh:startCh+numOfCh-1) + offset;

 for i = startCh:(startCh+numOfCh-1)
     channelPlot(:,i) = channel(:,i)+i*5000;
 end

%L = length(tmp);                      % Length of signal
time = (1:1e6)*T;
figure(1), plot(time,channelPlot);
  

for i = startCh:(startCh+numOfCh-1)
    ichannel(:,i) = ifft(abs(fft(channel(:,i))));
end
%figure(11), plot(ichannel);
hold on
%plot(channel);
hold off

%fft
figure(2);
for i = startCh:(startCh+numOfCh-1)
    tmp = fft(abs(channel(:,i)));
    fourier(:,i) = abs(tmp);
    subplot(2,2,i);
    plot(abs(tmp(1:(length(tmp)/2)+1))),title('Fourier:input') %Fs/2*linspace(0,1,(length(tmp)/2)+1),
    xlabel('Normalized frequency')
    ylabel('Magnitude')
end

%fft cleaning
figure(3);
for i = startCh:(startCh+numOfCh-1)
    tmp = cleanHarmonics(fourier(:,i));
    %tmp = cleanHarmonicsNeg(fourier(:,i));
    cleanFourier(:,i) = tmp;
    cleanChannel(:,i) = ifft(tmp);
    subplot(2,2,i);
    plot(abs(tmp(1:(length(tmp)/2)+1))),title('Fourier:cleaned') %Fs/2*linspace(0,1,(length(tmp)/2)+1),
    xlabel('Normalized frequency')
    ylabel('Magnitude')
end

%fft cleaning negatives
figure(33);
for i = startCh:(startCh+numOfCh-1)
    %tmp = cleanHarmonics(fourier(:,i));
    tmp = cleanHarmonicsNeg(fourier(:,i));
    cleanFourierNeg(:,i) = tmp;
    cleanChannel(:,i) = ifft(tmp); %CLEANING NEGATIVES
    subplot(2,2,i);
    plot(abs(tmp(1:(length(tmp)/2)+1))),title('Fourier:cleaned with negative') %Fs/2*linspace(0,1,(length(tmp)/2)+1),
    xlabel('Normalized frequency')
    ylabel('Magnitude')
end

 for i = startCh:(startCh+numOfCh-1)
     channelPlot(:,i) = ifft(cleanFourierNeg(:,i))+i*5000;
 end
figure(11), plot(time,channelPlot),title('cleaned inversed fourier channels');


%plot both;
figure(4);
for i = startCh:(startCh+numOfCh-1)
    tmp = fourier(:,i)+1000;
    subplot(2,2,i);
    plot(abs(tmp(1:(length(tmp)/2)+1))),title('Fourier:cleaned neg vs uncleaned') %Fs/2*linspace(0,1,(length(tmp)/2)+1),
    xlabel('Normalized frequency')
    ylabel('Magnitude')
    hold on;
    tmp = cleanFourierNeg(:,i);
    plot(abs(tmp(1:(length(tmp)/2)+1))),title('Fourier:cleaned neg vs uncleaned');
    hold off;
end

%plot both;
figure(44);
for i = startCh:(startCh+numOfCh-1)
    tmp = cleanFourier(:,i)+1000;
    subplot(2,2,i);
    plot(abs(tmp(1:(length(tmp)/2)+1))),title('Fourier:cleaned neg vs clean') %Fs/2*linspace(0,1,(length(tmp)/2)+1),
    xlabel('Normalized frequency')
    ylabel('Magnitude')
    hold on;
    tmp = cleanFourierNeg(:,i);
    plot(abs(tmp(1:(length(tmp)/2)+1))),title('Fourier:cleaned neg vs clean');
    hold off;
end

%running the high pass filter
for i = startCh:(startCh+numOfCh-1)
    filteredCleanChannel(:,i) = highPassFilter(cleanChannel(:,i), Fs);
end

for i = startCh:(startCh+numOfCh-1)
    filteredChannel(:,i) = highPassFilter(channel(:,i), Fs);
end

figure(5);
for i = startCh:(startCh+numOfCh-1)
    tmp = filteredCleanChannel(:,i);
    subplot(2,2,i);
    %plot(Fs/2*linspace(0,1,(length(tmp)/2)+1),tmp(1:(length(tmp)/2)+1)),title('filtered and cleaned channels') %Fs/2*linspace(0,1,(length(tmp)/2)+1),
    xlabel('Normalized frequency')
    ylabel('Magnitude')
end

% PLOT FILTERED CHANNELS
 for i = startCh:(startCh+numOfCh-1)
     channelPlot(:,i) = filteredCleanChannel(:,i)+i*5000;
 end
figure(55), plot(time,channelPlot),title('clean and filtered channels');

figure(6), plot(filteredCleanChannel);
hold on
plot(channel), title('filtered vs regular');
hold off

%figure(2), plot(abs(fourier(1:(length(fourier)/2)+1,:)));

%plot(fft(channel));p(end

%filtered uncleaned channel
figure(7);
for i = startCh:(startCh+numOfCh-1)
    tmp = filteredChannel(:,i);
    subplot(2,2,i);
    %plot(Fs/2*linspace(0,1,(length(tmp)/2)+1),abs(tmp(1:(length(tmp)/2)+1))),title('filtered') %Fs/2*linspace(0,1,(length(tmp)/2)+1),
    xlabel('Normalized frequency')
    ylabel('Magnitude')
end

% PLOT FILTERED CHANNELS
 for i = startCh:(startCh+numOfCh-1)
     channelPlot(:,i) = filteredChannel(:,i)+i*5000;
 end
figure(77), plot(time,channelPlot),title('filtered channels');


end


function tmp=cleanHarmonicsNeg(input)
offset = 2*1e4;
tmp=abs(input);%input(1:length(input)/2+1);%(offset:end);  %REMOVE!
out=input;
period = 100;
%mean100 = mean(out(1:100));
sumPer = sum(tmp(offset+1:offset+period)); 
for i = offset+period+1:length(tmp) %2005050
    if tmp(i)>3*sumPer/period       %20000
        tmp(i)=sumPer/period;
        out(i)= sumPer/period;
        if out(i)<0
            out(i)= -1*sumPer/period;
        end
    end
    sumPer=sumPer-tmp(i-period)+tmp(i);
end

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
% Butterworth Bandpass filter designed using FDESIGN.BANDPASS.
%for details
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

function plotting

figure(1);
plot(tmp,title('input'));
%https://www.mathworks.com/matlabcentral/answers/22119-fft-of-100-data-points
%{ 
L=length(x);
dt=3.5*10^(-6); %time.
fs=1/dt;
t=(0:1:L-1)*dt; %If you'd like to plot the time signal: plot(t,x)
out=fft(x,L)/L;
%}
%FOURIER
%normalizing and plotting
%plot(abs(fft(tmp)));
L = length(tmp);                      % Length of signal
time = (0:L-1)*T;                % Time vector
out=fft(tmp)/L;
figure(2), plot(Fs/2*linspace(0,1,(length(out)/2)+1),abs(out(1:(length(out)/2)+1))),title('Fourier: One sided Spectrum')
xlabel('Normalized frequency')
ylabel('Magnitude')
%figure(5), plot(abs(fftshift(fft(tmp))));



out=fft(filteredData)/L;
figure(3), plot(Fs/2*linspace(0,1,(length(out)/2)+1),abs(out(1:(length(out)/2)+1))),title('Filtered Fourier')
xlabel('Normalized frequency')
ylabel('Magnitude')

%PLOTTING
filtered = filtfilt(Hd.sosMatrix,Hd.ScaleValues,tmp);
figure(4);
plot(time,tmp-(mean(tmp)-mean(filtered)));
hold on
plot(time,filtered)
xlabel('Time (s)')
ylabel('Amplitude')
legend('Original Signal','Filtered Data');


plot(fail);
end



