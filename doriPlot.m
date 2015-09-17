clear all; close all; clc;
%{
2015-01-01_15-45-34_dreadd_rat_ref_animalground_2100depth     40 million samples
2015-01-29_dreadd_ref17_with_commutator_thresh50              15 million samples
2015-01-06_130_rat130_arena1                                  50 million samples
2015-05-18_messi_before_injection_threshold40_ref4            60 million samples
%}
sdir = {'2015-01-01_15-45-34_dreadd_rat_ref_animalground_2100depth'; 
         '2015-01-29_dreadd_ref17_with_commutator_thresh50'; 
         '2015-01-06_130_rat130_arena1' ;
         '2015-05-18_messi_before_injection_threshold40_ref4'};
file = char(sdir(1));
hdir = 'C:\\Users\\alm\\Desktop\\dori\\raw_data';
ADBitVolts = 0.000000015624999960550667;;
ch = 3;
depth = load(sprintf('%s\\%s\\ch%d.csv',hdir,file,ch));
depth = depth*ADBitVolts; % to volts
depth = depth*1e6; % to micro volts
tmp = depth(1:3e5);
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
Fs = 32000;                   % Sampling frequency
T = 1/Fs;                     % Sample time
L = length(tmp);                      % Length of signal
time = (0:L-1)*T;                % Time vector
out=fft(tmp)/L;
figure(2), plot(Fs/2*linspace(0,1,(length(out)/2)+1),abs(out(1:(length(out)/2)+1))),title('Fourier: One sided Spectrum')
xlabel('Normalized frequency')
ylabel('Magnitude')
%figure(5), plot(abs(fftshift(fft(tmp))));

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
%toDecibal = 20*log(10);         
Hd = design(fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs),'butter');
out=fft(filtfilt(Hd.sosMatrix,Hd.ScaleValues,tmp))/L;
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



%figure(5), plot(time,filtered);





plot(fail);

part = depth(1:L);

figure(1)
plot(Fs*t,part)
title('Original Signal')
xlabel('time')

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(part,NFFT)/L;
plot(Y);

% Plot single-sided amplitude spectrum.
figure(2)
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of signal')
xlabel('Frequency [Hz]')
ylabel('Amplitude [\mu V]')

% High Pass Filter
Fstop = 550;
Fpass = 600;
Astop = 65;
Apass = 0.5;

d = designfilt('highpassfir','StopbandFrequency',Fstop, ...
  'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
  'PassbandRipple',Apass,'SampleRate',Fs,'DesignMethod','equiripple');

% fvtool(d)

Y = fftfilt(d,part,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure(3)
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of signal')
xlabel('Frequency [Hz]')
ylabel('Amplitude [\mu V]')

Y_ymp = Y;
% Y_tmp(1:599) = 0;
Z = ifwft(Y_tmp);

figure(4)
plot(Fs*t,Z(1:L))
title('droped Signal')
xlabel('time')