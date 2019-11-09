clc
clear all
close all

A = xlsread('Damping.xlsx')

figure(1)
plot(A(:,1),A(:,2))

t = A(:,1);
x = A(:,2)
% 
% The signal is sampled at 250 MHz
FS = 250e6 ; % Sampling frequency
% According to Nyquist theorm we can only have frequencies upto 12.5 MHz
%   in this signal. The function butter take normalised frequency.
FN = FS/2; % Nyquist frequency

figure(1)
plot (t,x) ; grid
title ('Time Plot of Unfiltered Signal') 
xlabel('Time (sec)')
ylabel('Magnitude')

% Plot first half of DFT (normalised frequency)
X_mags = abs(fft(x))/FS;
num_bins = length(X_mags);

figure(2)
plot([0:1/(num_bins/2 -1):1]*FN, X_mags(1:(num_bins-1)/2))
title ('Frequency Plot of Unfiltered signal') 
xlabel('Frequency (Hz)')
ylabel('Magnitude')
format long

x_dfiltered=(ADC11_F.data)';
 
% Plot first half of DFT (normalised frequency)
X_mags = abs(fft(x_dfiltered))/FS;
num_bins = length(X_mags);

figure(4)
plot([0:1/(num_bins/2 -1):1]*FN, X_mags(1:(num_bins-1)/2))
title ('Frequency Plot of Filtered signal') 
xlabel('Frequency (Hz)')
ylabel('Magnitude')

figure(5)
plot (t,x_dfiltered) ; grid
title ('Time Plot of Filtered Signal') 
xlabel('Time (sec)')
ylabel('Magnitude')

% Design butterworth filter 
Wp=[171 213]./FN%Passband (normalised) corner frequencies (rad/sample) 
Ws=[172.8 211.2]./FN%Stopband (normalised) corner frequencies (rad/sample) 
Rp=3%Maximum Passband Ripple (dB)
Rs=-10*log10(1/5)%Minimum Attenuation of the Undesired Frequency range (dB)

[n,Wn] = buttord(Wp,Ws,Rp,Rs) %n: Order of the Filter    Wn: Normalized cutoff frequencies     
[b a] = butter(n, Wn, 'stop');
 
%plot the frequency response (normalised frequency)
H = freqz(b,a, floor(num_bins/2));

figure(3); 
plot([0:1/(num_bins/2 -1):1]*FN, abs(H),'r');grid
title ('Frequency Response of Difference Equation form Filter') 
xlabel('Frequency (Hz)')
ylabel('Magnitude') 

x_DEfiltered=filter(b, a, x);

% Plot first half of DFT (normalised frequency)
X_mags = abs(fft(x_DEfiltered))/FS;
num_bins = length(X_mags);

figure(6)
plot([0:1/(num_bins/2 -1):1]*FN, X_mags(1:(num_bins-1)/2))
title ('Frequency Plot of Difference Equation form Filtered signal') 
xlabel('Frequency (Hz)')
ylabel('Magnitude')

figure(7)
plot (t,x_DEfiltered) ; grid
title ('Time Plot of Difference Equation form Filtered Signal') 
xlabel('Time (sec)')
ylabel('Magnitude')


