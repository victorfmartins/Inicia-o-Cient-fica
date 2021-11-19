%% TerceiroMes - Phase-amplitude coupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear, clf, clc
srate = 1000; % in Hz
dt = 1/srate; % in s
t = dt:dt:20; % in s

%%% Toy Model Creation %%%
% set sine waves
nonModAmp = 10; % increase this to get less modulation (lower MI value) 
                % and lower artifacts in the ampfreq comodulogram
slow_modulation_wave = sin(2*pi*8*t);
fast_modulated_wave = sin(2*pi*80*t).*(0.2*(slow_modulation_wave+1)+nonModAmp*0.1);
% fast_modulated_wave = sin(2*pi*80*t).*slow_modulation_wave;

% fuse sine waves + noise to create signal
LFP = slow_modulation_wave + ...
      fast_modulated_wave  + ...
      1*randn(size(t));
  
%%% Feature Extraction %%%
% filter the signal for low freq for the phase
% modulating wave and for higth freq for the amp
% modulated wave
phase_freq = eegfilt(LFP,srate,5,10);
amp_freq = eegfilt(LFP,srate,30,100);

% all set to extract respective phases and
% amplitudes
phase = angle(hilbert(phase_freq));
ampenv = abs(hilbert(amp_freq));

%%% Capability demonstration %%%
% show we can select fractions of the wave with
% specific features

% each j is selecting a band with 20 specific
% phase points 
j=17; % 1 < j <= 18 bin de 20 pontos (360/20=18)
I = find(phase > deg2rad(-180+j*20) & ...
         phase < deg2rad(-160+j*20));

% get the mean point of the corresponding band in
% the amplitude envelop vector. This is the mean
% over all the points of all the many 20 point bands
% of the signal. That is, if the low freq wave
% oscilates at 5Hz, them in one second this will
% get the mean over 5*20 points
MeanAmp = mean(ampenv(I));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Figure 1 Tort et al. 2010 %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(111), clf
    plot(t,LFP,'k-'); hold on % plot the signal in black
    plot(t(I),LFP(I),'y.');
    plot(t,phase_freq-4,'b-','linew',2)  % plot the low pass filtered signal
    plot(t(I),phase_freq(I)-4,'y.','linew',2)
    plot(t,1*amp_freq-7,'r-','linew',0.5)% plot the high pass filtered signal
    plot(t,1*ampenv-7,'r-','linew',2) % plot its amplitude envelop 
    plot(t(I),1*ampenv(I)-7,'y.') % plot the selected band in the amplitude envelop
    plot(t,phase-12,'b.') % plot phase 
    plot(t(I),phase(I)-12,'y.'); hold off % plot the selected phase band
    ylim([-16 5])
    xlim([0 1])
    xlabel('Time (s)')
    ax = gca;
    ax.YTick = [-12 -8 -4 0];
    ax.YTickLabel = {'Theta Phase','Gamma','Theta','LFP'};
    title(['Mean \gamma Amplitude = ' num2str(MeanAmp)])

%% Modulation Index - MI
% the mean calculation in the last cell will be
% repeated over all 18 bins. This will create a 18
% points vector of the mean amplitude envelope of
% the hight freq wave. Its normalizatoin and
% entropy calculation is the MI metric calculation.

% For every bin take the mean of the amplitude envelope
clear MeanAmp
count = 0;
for phasebin = -180:20:160
    count = count+1;
    I = find(phase > deg2rad(phasebin) & ...
             phase < deg2rad(phasebin+20));
    MeanAmp(count) = mean(ampenv(I));
end

% Normalizing to get a deistribution of
% probability and entropy.
p = MeanAmp/sum(MeanAmp);

% computing entropy values:
H = -sum(p.*log(p));
Hmax = log(length(p));

% computing the MI metric
MI = (Hmax-H)/Hmax;

clf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Figure 1 D  Tort et al. 2010 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(121)
    bar(10:20:710,[MeanAmp MeanAmp])
    xlabel('\Theta Phase (Deg)')
    ylabel('Mean \gamma Amplitude')
    set(gca,'xtick',0:90:720)
    title(['Modulation Index = ' num2str(MI)])
    xlim([0 720])
subplot(122)
    bar(10:20:710,[p p])
    xlabel('\Theta Phase (Deg)')
    ylabel('Normalized mean \gamma Amplitude')
    set(gca,'xtick',0:90:720)
    title(['Modulation Index = ' num2str(MI)])
    xlim([0 720])

%%
% % Comodulograma

clear Comodulogram

phase_freq_vector = 0:2:20;
fp_bandwidth = 4;

amp_freq_vector = 20:5:200;
fa_bandwidth = 20;

count_phase = 0;
tic
for fp = phase_freq_vector
count_phase = count_phase+1
phase_freq = eegfilt(LFP,srate,fp,fp+fp_bandwidth);
phase = angle(hilbert(phase_freq));
count_amp = 0;
for fa = amp_freq_vector
    count_amp = count_amp+1;
    amp_freq = eegfilt(LFP,srate,fa,fa+fa_bandwidth);
    ampenv = abs(hilbert(amp_freq));
    % MI calculation
    count = 0;
    for phasebin = -180:20:160
        count = count+1;
        I = find(phase > deg2rad(phasebin) & ...
                 phase < deg2rad(phasebin+20));
        MeanAmp(count) = mean(ampenv(I));
    end
    p = MeanAmp/sum(MeanAmp);
    H = -sum(p.*log(p));
    Hmax = log(length(p));
    MI = (Hmax-H)/Hmax;
    Comodulogram(count_phase,count_amp) = MI;
    end
end
toc

subplot(111), clf
    contourf(phase_freq_vector+fp_bandwidth/2,...
               amp_freq_vector+fa_bandwidth/2,...
               Comodulogram',50,'linestyle','none')
    axis xy
    xlabel('Phase Frequency (Hz)')
    ylabel('Amplitude Frequency (Hz)')
    
        
%% Figure 2 Tort el al. 2010 - Visualizing Modulation

clear, clc, clf
srate = 1000; % in Hz
dt = 1/srate; % in s
t = dt:dt:20; % in s

ivector = [200 8 4.5 3];
for i = 0:3
%%% Toy Model Creation %%%
% set sine waves
nonModAmp = ivector(i+1); % increase this to get less modulation (lower MI value) 
                % and lower artifacts in the ampfreq comodulogram
slow_modulation_wave = sin(2*pi*8*t);
fast_modulated_wave = sin(2*pi*80*t).*(0.2*(slow_modulation_wave+1)+nonModAmp*0.1);
% fast_modulated_wave = sin(2*pi*80*t).*slow_modulation_wave;

% fuse sine waves + noise to create signal
LFP = slow_modulation_wave + ...
      fast_modulated_wave  + ...
      0.1*randn(size(t));
  
%%% Feature Extraction %%%
% filter the signal for low freq for the phase
% modulating wave and for higth freq for the amp
% modulated wave
phase_freq = eegfilt(LFP,srate,5,10);
amp_freq = eegfilt(LFP,srate,30,100);

% all set to extract respective phases and
% amplitudes
phase = angle(hilbert(phase_freq));
ampenv = abs(hilbert(amp_freq));

subplot(4,4,[1 2]+(4*i)),
    plot(t,LFP,'k-'); hold on % plot the signal in black
    m = max(LFP);
    plot(t,1*amp_freq-2.2*m,'r-','linew',0.5)% plot the high pass filtered signal
    plot(t,1*ampenv-2.2*m,'r-','linew',2); hold off % plot its amplitude envelop
    xlim([0 1])
    ylabel(['Case ' num2str(i+1)])
    if (i == 3)
        xlabel('Time (s)')
    end
    

clear MeanAmp
count = 0;
for phasebin = -180:20:160
    count = count+1;
    I = find(phase > deg2rad(phasebin) & ...
             phase < deg2rad(phasebin+20));
    MeanAmp(count) = mean(ampenv(I));
end

% Normalizing to get a deistribution of
% probability and entropy.
p = MeanAmp/sum(MeanAmp);

% computing entropy values:
H = -sum(p.*log(p));
Hmax = log(length(p));

% computing the MI metric
MI(i+1) = (Hmax-H)/Hmax;

subplot(4,4,3+(4*i))
    bar(10:20:710,[p p])
    set(gca,'xtick',0:360:720)
    xlim([0 720])
    ylim([0 0.1])
    ylabel('Amplitude')
    if (i == 3)
        xlabel('Phase (Deg)')
    end
end

subplot(4,4,[4 8 12 16])
    bar(MI)
    
%% Effect of epoch length and noise - Figure 3A (data gathering)

clear, clf, clc
srate = 1000; % in Hz
dt = 1/srate; % in s

tic
MI1 = zeros(6,100,100);
noisevector = [0 0.08 0.16 0.24 0.32 0.4];
for kk = 1:6 % mudando a intensidade do ruindo
kk
for jj = 1:100 % mudando o tempo maximo
jj
t = dt:dt:jj; % in s

%%% Toy Model Creation %%%
% set sine waves
nonModAmp = 10; %% should be 13
slow_modulation_wave = sin(2*pi*7.5*t);
fast_modulated_wave = sin(2*pi*45*t).*(0.2*(slow_modulation_wave+1)+nonModAmp*0.1);

MI = zeros(1,100);
for ii = 1:100 % repetindo o calculo 100 para gerar o CV
% fuse sine waves + noise to create signal
LFP = slow_modulation_wave + ...
      fast_modulated_wave  + ...
      noisevector(kk)*randn(size(t));

%%% Feature Extraction %%%
% filter the signal for low freq for the phase
% modulating wave and for higth freq for the amp
% modulated wave
phase_freq = eegfilt(LFP,srate,5,10,0,1);
amp_freq = eegfilt(LFP,srate,30,60,0,1);

% all set to extract respective phases and
% amplitudes
phase = angle(hilbert(phase_freq));
ampenv = abs(hilbert(amp_freq));

clear MeanAmp
count = 0;
for phasebin = -180:20:160
    count = count+1;
    I = find(phase > deg2rad(phasebin) & ...
             phase < deg2rad(phasebin+20));
    MeanAmp(count) = mean(ampenv(I));
end

% Normalizing to get a deistribution of
% probability and entropy.
p = MeanAmp/sum(MeanAmp);

% computing entropy values:
H = -sum(p.*log(p));
Hmax = log(length(p));

% computing the MI metric
MI(ii) = (Hmax-H)/Hmax;
MI1(kk,jj,ii) = (Hmax-H)/Hmax;
end
end
end
toc
save('C:\Users\VAITA\OneDrive\Avançado\Scientific-Research\PACdata\bigA_MI.mat', 'MI1');
% Elapsed time is  1.0258e+03 seconds.

    
%% Figure 3A - ploting

clear, clc, clf
% load('MIa.mat')
load('PACdata\bigA_MI.mat')
%  kk: noise;
%  jj: length;
%  ii: trial;

for i = 1:10
if (i<4 | (i>5 & i<9))
P1 = subplot(2,5,mod(i,5)+5*(i>5));
    mean_MI1 = mean(squeeze(MI1(mod(i,5)+3*(i>5),:,:))');
    std_MI1 = std(squeeze(MI1(mod(i,5)+3*(i>5),:,:))');
    color = [.9 .9 .9]-[.1 .1 .1]*(i-2*(i>5));
    plot((1:100), mean_MI1,'k','linew',2,'color',color), hold on
    plot((1:100), mean_MI1-std_MI1,'k--','linew',2,'color',color)
    plot((1:100), mean_MI1+std_MI1,'k--','linew',2,'color',color), hold off
    ylim([0 0.002]);
end
P1 = subplot(2,5,7);
xlabel('Epoch Length (s)')
ylabel('Modulation Index')
set( get(P1,'YLabel'), 'Position', [-180 0.0032] );
end
P2 = subplot(2,5,[4 5 9 10]);
    mean_MI1 = mean(MI1(:,:,:),3);
    std_MI1 = std(MI1(:,:,:),0,3);
    color = [(.9-.1*(1:6))' (.9-.1*(1:6))' (.9-.1*(1:6))'];
    hL = plot([(1:100)],(std_MI1./mean_MI1)','linew',2);
    for i=1:length(hL),set(hL(i),'color',color(i,:));end
%     set(hL,{'color'},mat2cell(color,length(color),1))
    xlabel('Epoch Length (s)')
    ylabel('Coefficient of Variation')
    noise = [0 0.08 0.16 0.24 0.32 0.4];
    legend('\sigma = 0','\sigma = 0.08','\sigma = 0.16','\sigma = 0.24','\sigma = 0.32','\sigma = 0.40')
    set( get(P2,'YLabel'), 'Position', [110 0.0700] );
    ylim([0 0.6])


%% plot first epoch lenght point to reach 10%CV as a function of noise
% run simulation until treshold only!!!!!!!!!!!

clear, clf, clc
srate = 1000; % in Hz
dt = 1/srate; % in s

tic
CVlim = zeros(1,75);
MI1 = zeros(75,200,100);
noisevector = 0:0.02:1.5-.02;
kk = 1;
while (kk < 75) % mudando a intensidade do ruindo

kk
if (kk > 1 & CVlim(kk-1) - 2 > 0)
    jj = CVlim(kk-1) - 2;
elseif (kk > 50)
    jj = CVlim(kk-1) + 2;
else
    jj = 1;
end

while (jj < 200) % mudando o tempo maximo
% if (jj > 1)
%     jj = CVlim(kk);
% end
jj
t = dt:dt:jj; % in s

%%% Toy Model Creation %%%
% set sine waves
nonModAmp = 10; %% should be 13
slow_modulation_wave = sin(2*pi*7.5*t);
fast_modulated_wave = sin(2*pi*45*t).*(0.2*(slow_modulation_wave+1)+nonModAmp*0.1);

MI = zeros(1,100);
for ii = 1:100 % repetindo o calculo 100 para gerar o CV
% fuse sine waves + noise to create signal
LFP = slow_modulation_wave + ...
      fast_modulated_wave  + ...
      noisevector(kk)*randn(size(t));

%%% Feature Extraction %%%
% filter the signal for low freq for the phase
% modulating wave and for higth freq for the amp
% modulated wave
phase_freq = eegfilt(LFP,srate,5,10,0,1);
amp_freq = eegfilt(LFP,srate,30,60,0,1);

% all set to extract respective phases and
% amplitudes
phase = angle(hilbert(phase_freq));
ampenv = abs(hilbert(amp_freq));

clear MeanAmp
count = 0;
for phasebin = -180:20:160
    count = count+1;
    I = find(phase > deg2rad(phasebin) & ...
             phase < deg2rad(phasebin+20));
    MeanAmp(count) = mean(ampenv(I));
end

% Normalizing to get a deistribution of
% probability and entropy.
p = MeanAmp/sum(MeanAmp);

% computing entropy values:
H = -sum(p.*log(p));
Hmax = log(length(p));

% computing the MI metric
MI(ii) = (Hmax-H)/Hmax;
MI1(kk,jj,ii) = (Hmax-H)/Hmax;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_MI1 = mean(MI1(kk,jj,:),3);
std_MI1 = std(MI1(kk,jj,:),0,3);
CV = (std_MI1./mean_MI1);
if (CV <= 0.1)
    CVlim(kk) = jj;
    break;
end
jj = jj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
kk = kk + 1;
end
toc
save('C:\Users\VAITA\OneDrive\Avançado\Scientific-Research\PACdata\CVlim', 'CVlim');


%% Plot

clear, clc, clf
load('PACdata\CVlim.mat')
noisevector = 0:0.02:1.5-.02;
plot(noisevector,CVlim)
    xlabel('Noise \sigma')
    ylabel('Epoch Length When Crossing CV < 10% Treshold (s)')
    title('Size of Epoch for CV<10%')
%     xlim([0 1.5])

%% Effect of epoch length and coupling strength - Figure 3B

clear, clf, clc
srate = 1000; % in Hz
dt = 1/srate; % in ms

tic
MI2 = zeros(6,100,100);
for kk = 1:6 % mudando a intensidade do acoplamento
kk
for jj = 1:100 % mudando o tempo maximo
t = dt:dt:jj; % in s

%%% Toy Model Creation %%%
% set sine waves
nonModAmp = [27 22 18 15 13 12]; % quanto maior menor é coupling
slow_modulation_wave = sin(2*pi*7.5*t);
fast_modulated_wave = sin(2*pi*45*t).*(0.2*(slow_modulation_wave+1)+nonModAmp(kk)*0.1);

MI = zeros(1,100);
for ii = 1:100 % repetindo o calculo 100 para gerar o CV
% fuse sine waves + noise to create signal
LFP = slow_modulation_wave + ...
      fast_modulated_wave  + ...
      0.24*randn(size(t));

%%% Feature Extraction %%%
% filter the signal for low freq for the phase
% modulating wave and for higth freq for the amp
% modulated wave
phase_freq = eegfilt(LFP,srate,5,10,0,1);
amp_freq = eegfilt(LFP,srate,30,60,0,1);

% all set to extract respective phases and
% amplitudes
phase = angle(hilbert(phase_freq));
ampenv = abs(hilbert(amp_freq));

clear MeanAmp
count = 0;
for phasebin = -180:20:160
    count = count+1;
    I = find(phase > deg2rad(phasebin) & ...
             phase < deg2rad(phasebin+20));
    MeanAmp(count) = mean(ampenv(I));
end

% Normalizing to get a deistribution of
% probability and entropy.
p = MeanAmp/sum(MeanAmp);

% computing entropy values:
H = -sum(p.*log(p));
Hmax = log(length(p));

% computing the MI metric
MI(ii) = (Hmax-H)/Hmax;
MI2(kk,jj,ii) = (Hmax-H)/Hmax;
end
end
end
toc
save('C:\Users\VAITA\OneDrive\Avançado\Scientific-Research\PACdata\bigB_MI.mat', 'MI2');
% Elapsed time is 841.489517 seconds.
 

%% Figure 3B - ploting

clear, clc, clf
load('PACdata\bigB_MI.mat')
%  kk: modAmp;
%  jj: length;
%  ii: trial;

for i = 1:10
if (i<4 | (i>5 & i<9))
P1 = subplot(2,5,mod(i,5)+5*(i>5));
    mean_MI2 = mean(squeeze(MI2(mod(i,5)+3*(i>5),:,:))');
    std_MI2 = std(squeeze(MI2(mod(i,5)+3*(i>5),:,:))');
    color = [.9 .9 .9]-[.1 .1 .1]*(i-2*(i>5));
    plot((1:100), mean_MI2,'k','linew',2,'color',color), hold on
        plot((1:100), mean_MI2-std_MI2,'k--','linew',2,'color',color)
    plot((1:100), mean_MI2+std_MI2,'k--','linew',2,'color',color), hold off
    ylim([0 0.001]);
end
P1 = subplot(2,5,7);
xlabel('Epoch Length (s)')
ylabel('Modulation Index')
set( get(P1,'YLabel'), 'Position', [-180 0.0032] );
end
P2 = subplot(2,5,[4 5 9 10]);
    mean_MI2 = mean(MI2(:,:,:),3);
    std_MI2 = std(MI2(:,:,:),0,3);
    color = [(.9-.1*(1:6))' (.9-.1*(1:6))' (.9-.1*(1:6))'];
    hL = plot([(1:100)],(std_MI2./mean_MI2)','linew',2);
    for i=1:length(hL),set(hL(i),'color',color(i,:));end
%     ylim([0 .8])
    xlabel('Epoch Length (s)')
    ylabel('Coefficient of Variation')
    legend('Case 1','Case 2','Case 3','Case 4','Case 5','Case 6')
    set(get(P2,'YLabel'),'Position',[110 0.275]);


%% Rate of response to Non Modulation Amplitude

clear, clf, clc
srate = 1000; % in Hz
dt = 1/srate; % in ms

tic
MI2 = zeros(50,1,100);
for kk = 1:50 % mudando a intensidade do acoplamento
kk
jj = 100; % mudando o tempo maximo
t = dt:dt:jj; % in s

%%% Toy Model Creation %%%
% set sine waves
nonModAmp = kk; % quanto maior menor é coupling
slow_modulation_wave = sin(2*pi*7.5*t);
fast_modulated_wave = sin(2*pi*45*t).*(0.2*(slow_modulation_wave+1)+nonModAmp*0.1);

MI = zeros(1,100);
for ii = 1:100 % repetindo o calculo 100 para gerar o CV
% fuse sine waves + noise to create signal
LFP = slow_modulation_wave + ...
      fast_modulated_wave  + ...
      0.24*randn(size(t));

%%% Feature Extraction %%%
% filter the signal for low freq for the phase
% modulating wave and for higth freq for the amp
% modulated wave
phase_freq = eegfilt(LFP,srate,5,10,0,1);
amp_freq = eegfilt(LFP,srate,30,60,0,1);

% all set to extract respective phases and
% amplitudes
phase = angle(hilbert(phase_freq));
ampenv = abs(hilbert(amp_freq));

clear MeanAmp
count = 0;
for phasebin = -180:20:160
    count = count+1;
    I = find(phase > deg2rad(phasebin) & ...
             phase < deg2rad(phasebin+20));
    MeanAmp(count) = mean(ampenv(I));
end

% Normalizing to get a deistribution of
% probability and entropy.
p = MeanAmp/sum(MeanAmp);

% computing entropy values:
H = -sum(p.*log(p));
Hmax = log(length(p));

% computing the MI metric
MI(ii) = (Hmax-H)/Hmax;
MI2(kk,jj,ii) = (Hmax-H)/Hmax;
end
end
toc
save('C:\Users\VAITA\OneDrive\Avançado\Scientific-Research\PACdata\nonModAmpEff.mat', 'MI2');
% Elapsed time is 289.070380 seconds.

%% Plot rate of response to Non Modulation Amplitude
% É linear. 
% Quanto maior o nonModAmp maior o CV.
% Quanto maior a modulação menor é o CV.

clear, clc, clf
load('PACdata\nonModAmpEff.mat')
%  kk: modAmp; 1:50
%  jj: length; 100
%  ii: trial;  1:100

subplot(111)
    mean_MImodAmp = mean(MImodAmp(:,:,:),3);
    std_MImodAmp = std(MImodAmp(:,:,:),0,3);
%     color = [(.9-.1*(1:6))' (.9-.1*(1:6))' (.9-.1*(1:6))'];
    hL = plot([(1:50)],(std_MImodAmp./mean_MImodAmp)','linew',2);
%     for i=1:length(hL),set(hL(i),'color',color(i,:));end
    xlabel('Non Modulation Amplitude')
    ylabel('Coefficient of Variation')


%%    
vv = [1 2 3 6 7 8 11 12 13 16 17 18];
count = 0;
for ww = 1:10 % change subplot
if (ww in vv)
count = count + 1;
subplot(5,4,ww)
    plot((1:25)*4,CVAll(count,:),'linew',2)
end
end

% Elapsed time is  1.0258e+03 seconds.
subplot(5,4,[4 5 9 10])
    plot((1:25)*4,CVAll,'linew',2)
    xlabel('Epoch Length (s)')
    ylabel('Coefficient of Variation')
    legend('Noise = 1','Noise = 2','Noise = 3','Noise = 4','Noise = 5')


































