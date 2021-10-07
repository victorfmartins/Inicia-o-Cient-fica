%% Code for Report 1 Figures

%% Figure 2 - Exemple of neuronal time series and theta band
clear, clc, clf
load('EE.210.mat')
Tmax = length(data)/sampleRate;
dt = 1/sampleRate;
t = dt:dt:Tmax;

figure(2), clf
subplot(2,1,1)
    plot(t-30, data(1,:), 'k')
    xlim([0 5])
    xlabel('Time (s)')
    ylabel('Voltage (mV)')
    title('Raw signal')
    filtrado = eegfilt(data(1,:),sampleRate,2,6);
    subplot(2,1,2)
    plot(t-30, filtrado(:), 'k')
    xlim([0 5])
    xlabel('Time (s)')
    ylabel('Voltage (mV)')
    title('Bandpass Filtered (8-12hz)')

%% Figure 3 - Coherence Spectrum

clear, clf, clc,
srate = 1000;
dt = 1/srate;
% changing the sample size changes the coherence
% as it changes the amount of vectors I am
% computing the average over
t = dt:dt:25; 

% the value randn add makes no difference as it is
% not changing over time 
phi = -deg2rad(90)+0.5*randn;

% if we remove the randn from X and Y we will have
% coherence 1 in all frequencies this is what is
% happening in the program  
X = 3*sin(2*pi*4*t)+5*randn(size(t));
Y = 1*sin(2*pi*4*t+phi)+3*randn(size(t));

% increasing window size generates more precision
% in frequencies but generates more noise in
% others.  
windowlength = 5*srate;
overlap = 0;

% freqvector = 0:0.1:20;
% [Cxy F] = mscohere(X,Y,windowlength,overlap,freqvector,srate)
% plot(F,Cxy) 
% xlabel('Frequency (Hz)')
% ylabel('Coherence')

nfft = 2^16;

figure(3), clf
subplot(4,1,1)
    plot(t,X), hold on
    plot(t,Y - 20), hold off
    xlim([0 2])
    ylim([-40 20])
    xlabel('Time (s)')
    title('Raw signal')
    legend('Signal A', 'Signal B')
subplot(4,1,2)
    Xfilt = eegfilt(X,srate,2,6);
    Yfilt = eegfilt(Y,srate,2,6);
    plot(t,Xfilt, 'linew', 2), hold on
    plot(t,Yfilt, 'linew', 2), hold off
    xlim([0 2])
    xlabel('Time (s)')
    ylabel('Voltage (mV)')
    title('Bandpass Filtered (2-6hz)')
subplot(4,1,3)
    Xfilt = eegfilt(X,srate,6,10);
    Yfilt = eegfilt(Y,srate,6,10);
    plot(t,Xfilt, 'linew', 2), hold on
    plot(t,Yfilt, 'linew', 2), hold off
    xlim([0 2])
    xlabel('Time (s)')
    ylabel('Voltage (mV)')
    title('Bandpass Filtered (6-10hz)')
subplot(4,1,4)
    [Cxy, F] = mscohere(X,Y,windowlength,overlap,nfft,srate);
    plot(F,Cxy.^2, 'k', 'linew', 2)
    xlabel('Frequency (Hz)')
    ylabel('Coherence')
    title('Coherence Spectrum')
xlim([0 10])














