%clear all; close all; clc;
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
ADBitVolts = 0.000000015624999960550667;
depth = load(sprintf('%s\\%s\\ch%d.csv',hdir,file,1));
depth = depth*ADBitVolts;
figure(2);
plot(ch1l);

plot(t);
dft = DFT_Fun(t);

tnn = fft(t);