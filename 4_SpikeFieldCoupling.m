%% QuartoMes - Spike-Phase coupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Toy Model Creation %%%
% setup
clear, clc, clf
f=8;
srate=1000;
dt=1/srate;
t=dt:dt:10;

% set sine waves
slow_modulation_wave = sin(2*pi*8*t);
fast_modulated_wave = sin(2*pi*80*t).*slow_modulation_wave;

% fuse sine waves + noise to create signal
LFP = slow_modulation_wave + ...
      fast_modulated_wave  + ...
      1*randn(size(t));

% filter the signal for low freq for the phase
% modulating wave and for higth freq for the amp
% modulated wave
thetha = eegfilt(LFP,srate,5,10);
% gamma = eegfilt(LFP,srate,30,100);

% set index of periodic spike signal with
% consistent relation with the sine wave in the 8
% Hz LFP
th_phase = angle(hilbert(thetha));
Idxmodulado = find(th_phase > deg2rad(260-90) & ...
                   th_phase < deg2rad(270-90));
Idxmodulado =  Idxmodulado(randi(length(Idxmodulado),[1,100]));
% Plus some random spike indices
Idxruido = randi(length(LFP),[1,200]);
Idx = [Idxmodulado, Idxruido];

% assign 1 to every point were there is a spike
spkind = zeros(size(LFP));
spkind(Idx)=1;

subplot(311)
    plot(t,LFP); hold on
    plot(t(Idxmodulado),LFP(Idxmodulado),'ko');
    plot(t(Idxruido),LFP(Idxruido),'ro'); hold off
    xlim([0 2])
    xlabel('Time (s)')
    ylabel('mv')
    title('Spike activity over LFP')
subplot(312)
    plot(t,thetha); hold on
    plot(t(Idxmodulado),thetha(Idxmodulado),'ko');
    plot(t(Idxruido),thetha(Idxruido),'ro'); hold off
    xlim([0 2])
    xlabel('Time (s)')
    ylabel('mv')
    title('Spike activity over Theta Band')
subplot(313)
    plot(t(Idxmodulado),1*ones(size(Idxmodulado)),'ro','markerfacecolor','w')
    hold on
    plot(t(Idxruido),2*ones(size(Idxruido)),'bo','markerfacecolor','w')
    hold off
    ax = gca;
    ax.YTick = [1 2];
    ax.YTickLabel = {'Neuron Modulado','Neuron Ruido'};    
    xlabel('Time (s)')
    ylim([0.5 2.5])
    xlim([0 2])
    title('Rastergram')
    
    
%%
% vector of the frequencies in which PLC and Kappa
% will be taken
freqvector = 1:1/2:50;

clear Kappaspectrum PLVspectrum MIspectrum
tic
count = 0;
for f = freqvector
    count = count+1

    % filter signal to each freq band and get its phase
    filtrado = eegfilt(LFP,srate,f,f+4);
    phase = angle(hilbert(filtrado));

    %compute PLV and Kappa values
    [mu,PLV,sigma,CI,kappa]=anglemean(phase(Idx));

    % store the PLV and Kappa values
    PLVspectrum(count) = PLV;
    Kappaspectrum(count) = kappa;

    % compute MI and store it
    [counts, phasebins] = hist(phase(Idx),deg2rad(-170:20:170));
    p = counts/sum(counts);
    MI = (log(length(phasebins))+sum(p(p>0).*log(p(p>0))))/log(length(phasebins));
    MIspectrum(count) = MI;
end
toc

subplot(221)
    plot(freqvector+2,PLVspectrum)
    xlabel('Freq (Hz)')
    ylabel('PLV','fontsize',20)
    title(['PLV_{max} = ' num2str(max(PLVspectrum))])

subplot(222)
    plot(freqvector+2,Kappaspectrum)
    xlabel('Freq (Hz)')
    ylabel('\kappa','fontsize',20)
    title(['\kappa_{max} = ' num2str(max(Kappaspectrum))])

% The proportional distance of maximu point with
% next local minimum is higger in MI than in PLV
% and kappa

% MI is the only metric that gets more precise if
% one increase the frequency resolution. 8.3% gain
% from precision of 1Hz to precision of .25Hz. No
% precision improvement after that
subplot(223)
    plot(freqvector+2,MIspectrum)
    xlabel('Freq (Hz)')
    ylabel('MI','fontsize',20)
    title(['MI_{max} = ' num2str(max(MIspectrum))])

% One can see that MI goes linear with PLV^2 and kappa^2 
subplot(224)
    plot(Kappaspectrum,MIspectrum,'ko'); hold on
    plot(PLVspectrum,MIspectrum,'ro'); hold off
    xlabel('\kappa (black) and PLV (red)')
    ylabel('MI','fontsize',20)


%% Testar para modulação em função de FR














































