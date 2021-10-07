%% README
% This script studies the fourier transform.
% The ft() is aplied to real and complex signals.
% The ft() is aplied to real data.

%% Fourier with complex signal.
% neuronal signals are not complex!

clear, clf, clc
srate = 100; % sampling rate in Hz

% setup: 
% frequency, phase-lag, and amplitude vectors
% they are used to create a composite sinewave
t = -1:1/srate:1; % in Seconds
frex = [1, 3, 4, 7, 8, 11, 14]; % in Hz
phaseLag = [pi, pi/6, pi/4, pi/5, pi/10, pi/3, 0]; %in radians
amplit = [2, 1, 1/3, 1/5, 1/2, 1/3, 1];

% allocates space to sine waves
sine_waves1 = zeros(length(frex), length(t));
sine_waves2 = zeros(length(frex), length(t));

% Create sines and cossines for real and imaginary
% part of signal
for fi = 1:length(frex)
    sine_waves1(fi,:) = amplit(fi)*sin(2*pi*frex(fi).*t+phaseLag(fi));
end

for fi = 1:length(frex)
    sine_waves2(fi,:) = amplit(fi)*cos(2*pi*frex(fi).*t+phaseLag(fi));
end

%make complex wave
sine_waves = sine_waves1+1i*sine_waves2;

% make composed wave
signal = sum(sine_waves);

% set fourier variables
N = length(signal);
nyquist = srate/2;
fourierTime = ((1:N)-1)/N;
fourierCoefs = zeros(size(signal));
F = linspace(0,nyquist,floor(N/2)+1);

% compute the fourier transform
for fi = 1:N
    fourierSine = exp(-1i*2*pi*(fi-1).*fourierTime);
    fourierCoefs(fi) = sum(signal.*fourierSine);
end

fourierCoefs = fourierCoefs/N;
matlab_fft = fft(signal)/N;

% flip fourier for ploting
% when the signal is complex the
% coefficients are backwards  
fourierCoefsBackyards = flip(fourierCoefs); 
matlab_fftBackyards = flip(matlab_fft);

figure(1), clf
subplot(221)
    plot(t, real(exp( -2*pi*1i*(10).*fourierTime )))
    xlabel('time (a.u.)'), ylabel('Amplitude')
    title('Sine wave')
    
subplot(222)
    plot(t, real(signal), 'b')
    hold on
    % why the imag part is different from the real?
    plot(t, imag(signal), 'b:')
    title('Data')

subplot(212)
% Why don't you have to multiply by 2? 
% Because the signal is complex and so the fft is not mirrored.
    plot(F, abs(fourierCoefsBackyards(1:length(F))),'*-') 
    hold on
    plot(F, abs(matlab_fftBackyards(1:length(F))))
    xlabel('Frequency (Hz)')
    ylabel('Power (\muV)')
    title('Power spectrum derived from discrete Fourier transform')


%% Phase Coherence with complex signal

% Setup
clear, clc, clf
srate = 500; % in Hz
f = 4; % in Hz
t = -1:1/srate:1; % in Seconds
phaseLag = pi/2+pi/6; % in radians
N = length(t);
nyquist = srate/2;
F = linspace(0,nyquist,floor(N/2)+1);

% Complex sine wave
X = exp(-1i*2*pi*f.*t); %the - makes the real one lead
sine_waveB1 = sin(2*pi*f.*t + phaseLag);
sine_waveB2 = cos(2*pi*f.*t + phaseLag);
Y = sine_waveB1+1i*sine_waveB2;

% Fourier
Fx = fft(X);
Fy = fft(Y);

% Cross-Spectrum
Fxy = mean(Fx.*conj(Fy));
Fxx = mean(Fx.*conj(Fx));
Fyy = mean(Fy.*conj(Fy));

% Normalized Cross-Spectrum
C = Fxy/(sqrt((Fxx.*Fyy)));

% Coherence
Coh = abs(C);

fourierABackyards = flip(Fx); 

% Ploting
subplot(3,1,1)
    plot(t, real(X), 'b')
    hold on
    plot(t, imag(X), 'b:')
    hold off 

subplot(3,1,2)
    plot(t, real(Y), 'k')
    hold on
    plot(t, imag(Y), 'k:')
    hold off 

subplot(3,1,3)
    plot(F, abs(fourierABackyards(1:length(F)))/N, 'b')
    xlabel('Frequency (Hz)')
    ylabel('Power (\muV)')
    title('Power spectrum derived from FFT')


%% Phase Coherence with real (aka normal) signal
% Por que o espectro de coerencia resultante é diferente?
% Como a fft() não sofre com o principio da incerteza?


% Setup
clear, clf, clc, close
srate = 100; f = 16; nyquist = srate/2; % in Hz /original with srate=1000
dt = 1/srate; Tmax = 4; t = 0:dt:Tmax; % in Seconds

N = length(t);
F = linspace(0,nyquist,floor(N/2)+1);

%%% Nolte with many trials %%%
tic
Nwin = 25;
Fx = zeros(Nwin,N);
Fy = zeros(Nwin,N);
for nwin = 1:Nwin
    rng(nwin);
    phaseLag = deg2rad(90)+0.3*randn;
    X = sin(2*pi*f.*t)+0.3*randn(size(t));
    Y = sin(2*pi*f.*t + phaseLag)+0.3*randn(size(t));
    Fx(nwin,:) = fft(X); % Fx is a matrix over epochs
    Fy(nwin,:) = fft(Y);
end

%Cross-Spectrum
Fxy = mean(Fx.*conj(Fy));
Fxx = mean(Fx.*conj(Fx));
Fyy = mean(Fy.*conj(Fy));

% Normalized Cross-Spectrum
C = Fxy./(sqrt((Fxx.*Fyy)));

% Coherence
Coh = abs(C);
toc % Elapsed time is 0.019512 seconds.

%{
% %%% Nolte modificado %%%
% % Fourier
% Fx = fft(X); % AmX.*exp(-1i.*(PhaseX))
% Fy = fft(Y); % AmpY.*exp(-1i.*(PhaseY))
% %Cross-Spectrum
% Fxy = Fx.*conj(Fy); % AmpX.*AmpY.*exp(-1i.*(PhaseX-PhaseY))
% Fxx = Fx.*conj(Fx); % AmpX.^2.*exp(-1i.*(PhaseX-Phasex)) = AmpX.^2
% Fyy = Fy.*conj(Fy); % AmpY.^2
% %normalization factor
% nfFxy = (Fxx.*Fyy).^(1/2); % (AmpX.^2.*AmpY.^2).^(1/2) = AmpX.*AmpY
% % Normalized Cross-Spectrum
% C = Fxy./nfFxy; % exp(-1i.*(PhaseX-PhaseY))
% % Coherence
% Coh = abs(C); % ones(1,length(C))
%}
%{
% %%% Cohen Fourier %%%
% N=length(t);
% fourierTime = ((1:N)-1)/N;
% fourierCoefs = zeros(size(X));
% frequencies = linspace(0,nyquist,floor(N/2)+1); 
% for fi=1:N
%     % create sine wave for this frequency
%     fourierSine = exp( -1i*2*pi*(fi-1).*fourierTime );
%     % compute dot product as sum of point-wise elements
%     fourierCoefs(fi) = sum(fourierSine .* X);
% end
% % scale Fourier coefficients to original scale
% fourierCoefs = fourierCoefs / N; %%% igual a Fx
%}
%{
% %%% Aula 6 modificada %%%
% for ff = 1:N
%     K = exp(-1i*2*pi*ff*t);
%     W = hamming(length(t))';
%     FX = mean((W.*X).*K); % AmpX*exp(-1i*PhaseX)
%     FY = mean((W.*Y).*K); % AmpY*exp(-1i*PhaseY)
%     nFX = FX/abs(FX); % exp(-1i*PhaseX)
%     nFY = FY/abs(FY); % exp(-1i*PhaseY)
%     nFXY = nFX*conj(nFY); % exp(-1i*(PhaseX-PhaseY))
%     Cxy = abs(mean(nFXY)); % 1
%     CxySpectrum(ff) = Cxy; % ones(1,length(C)) %%% igual Coh
% end
%}
%
%%% Aula 6 freq first than win %%%
tic
nvetor = 25;
Vall = zeros(nwin,N);
for nwin = 1:nvetor
    for ff = 0:N-1
        K = exp(-1i*2*pi*ff*t);
        rng(nwin);
        phaseLag = deg2rad(90)+0.3*randn;
        X = sin(2*pi*f*t)+0.3*randn(size(t));
        Y = sin(2*pi*f*t+phaseLag)+0.3*randn(size(t));
        W = hamming(length(t))';
        FX = mean((W.*X).*K); % AmpX*exp(-1i*PhaseX)
        FY = mean((W.*Y).*K); % AmpY*exp(-1i*PhaseY)
        nFX = FX/abs(FX); % exp(-1i*PhaseX)
        nFY = FY/abs(FY); % exp(-1i*PhaseY)
        nFXY = nFX*conj(nFY); % exp(-1i*(PhaseX-PhaseY)))
        Vall(nwin, ff+1) = nFXY;
    end
end
CxySpectrum = abs(mean(Vall));
toc % Elapsed time is 11.797561 seconds.
%
%{
%%% Aula 6 %%%
tic
clear CxySpectrum
for ff = 0:N-1
    K = exp(-1i*2*pi*ff*t);
    nvetor = 25;
    clear Vall
    for nwindow = 1:nvetor %fake signal de N*1seg.
        rng(nwindow);
        phaseLag = deg2rad(90)+0.3*randn;
        X = sin(2*pi*f*t)+0.3*randn(size(t));
        Y = sin(2*pi*f*t+phaseLag)+0.3*randn(size(t));
        W = hamming(length(t))';
        FX = mean((W.*X).*K); % AmpX*exp(-1i*PhaseX)
        FY = mean((W.*Y).*K); % AmpY*exp(-1i*PhaseY)
        nFX = FX/abs(FX); % exp(-1i*PhaseX)
        nFY = FY/abs(FY); % exp(-1i*PhaseY)
        nFXY = nFX*conj(nFY); % exp(-1i*(PhaseX-PhaseY)))
        Vall(nwindow) = nFXY;
    end
    Cxy = abs(mean(Vall));
    CxySpectrum(ff+1) = Cxy;
end
toc % Elapsed time is 11.797561 seconds.
%}

% Ploting
subplot(3,1,1)
    plot(t, X)
    hold on
    plot(t, Y)
    xlabel('Time (s)')
    ylabel('Voltagem (mV)')
    title('Signals')

subplot(3,1,2)
    plot(F, abs(Fx(1,1:length(F)))*2/N) 
    hold on
    plot(F, abs(Fy(1,1:length(F)))*2/N)
    xlabel('Frequency (Hz)')
    ylabel('Power (\muV)')
    title('Power spectrum derived from FFT')

subplot(3,1,3)
    plot(F,CxySpectrum(1:floor(N/2)+1).^2)
    hold on
    plot(F,Coh(1:floor(N/2)+1).^2)
    xlabel('Frequency (Hz)')
    ylabel('Coherence')
    legend('Cxy^2', 'Coh^2')
    
ICoh = find(Coh>0.8);
ICxy = find(CxySpectrum>0.8);


%% Phase Coherence with actual data

% Setup: Real data
clear, clf, clc
load('EE.210.mat');
srate = sampleRate;
window = 2; % in seconds
signalA = data(2,50000:50000+2*window*srate);
signalB = data(3,50000:50000+2*window*srate);
t = -window:1/srate:window; % in Seconds
N = length(t);
nyquist = srate/2;
F = linspace(0,nyquist,floor(N/2)+1);

% Setup do sinal adicional
freqA = 10;
freqB = 10;
phaseLag = pi/6;

% Sine_wave emitida pelas atividades neurais A e B nas respectivas regiões
% sob os eletrodos A e B. Frequencia de A, freqA = freqB
X = sin(2*pi*freqA.*t);
Y = sin(2*pi*freqB.*t + phaseLag);

% Junção do dado real e do dado produzido
amp = (mean(signalA) + mean(signalB))/2;
signalA = signalA + amp*X;
signalB = signalB + amp*Y;

% Fourier
Fx = fft(signalA);
Fy = fft(signalB);

%Cross-Spectrum
Fxy = mean(Fx.*conj(Fy));
Fxx = mean(Fx.*conj(Fx));
Fyy = mean(Fy.*conj(Fy));

% Normalized Cross-Spectrum
C = Fxy/(sqrt((Fxx.*Fyy)));

% Coherence
Coh = abs(C);

% Ploting
subplot(3,1,1)
    plot(t, signalA, 'b')

subplot(3,1,2)
    plot(t, signalB, 'b')

subplot(3,1,3)
    h = polar([0 angle(C)],[0 2*abs(C)],'k');
    set(h,'linewidth',5)
    title('Coherence', 'FontSize', 9);


%% O que está dando errado com a coerência??

% por que a parte imaginária está dando sempre zero quando o sinal é real?
% por que seu valor é igual ao cosseno de phaseLag?