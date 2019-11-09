% The signal is sampled at 1000 Hz
FS = 1000 ; % Sampling frequency
% According to Nyquist theorm we can only have frequencies upto 500 Hz
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

% Design butterworth filter 
Wn=[172.8 211.2]./FN
[b a] = butter(10, Wn, 'stop')
 
%plot the frequency response (normalised frequency)
H = freqz(b,a, floor(num_bins/2));

figure(3); 
plot([0:1/(num_bins/2 -1):1]*FN, abs(H),'r');grid
title ('Frequency Response of Designed Filter') 
xlabel('Frequency (Hz)')
ylabel('Magnitude') 

%filter the signal using the b and a coefficients obtained from
%the butter filter design function
x_filtered = filter(b,a,x);
 
% Plot first half of DFT (normalised frequency)
X_mags = abs(fft(x_filtered))/FS;
num_bins = length(X_mags);

figure(4)
plot([0:1/(num_bins/2 -1):1]*FN, X_mags(1:(num_bins-1)/2))
title ('Frequency Plot of Filtered signal') 
xlabel('Frequency (Hz)')
ylabel('Magnitude')

figure(5)
plot (t,x_filtered) ; grid
title ('Time Plot of Filtered Signal') 
xlabel('Time (sec)')
ylabel('Magnitude')
