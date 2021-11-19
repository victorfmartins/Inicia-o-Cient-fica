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
    
    
%% Spectrum of MI, PLV, and Kappa
% vector of the frequencies in which MI, PLV, and Kappa
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
    [counts, phasebinsC] = hist(phase(Idx),deg2rad(-170:20:170));
    p = counts/sum(counts);
    MI = (log(length(phasebinsC))+sum(p(p>0).*log(p(p>0))))/log(length(phasebinsC));
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


%% Testar para modulação em função de FR in given phase bin
% alterar a probabilidade de um spike
% cair em um determinado fase bin variando
% a quantidade total de spikes


%%% Toy Model Creation %%%
% setup
clear, clc, clf
f=8;
srate=1000;
dt=1/srate;
t=dt:dt:20;

% set sine waves
nonModAmp = 10;
slow_modulation_wave = sin(2*pi*7.5*t);
fast_modulated_wave = sin(2*pi*45*t).*(0.2*(slow_modulation_wave+1)+nonModAmp*0.1);

% fuse sine waves + noise to create signal
LFP = slow_modulation_wave + ...
      fast_modulated_wave  + ...
      0.24*randn(size(t));

% filter the signal for low freq for the phase
% modulating wave and for higth freq for the amp
% modulated wave
thetha = eegfilt(LFP,srate,5,10);
% gamma = eegfilt(LFP,srate,30,100);

% set index of periodic spike signal with
% consistent relation with the sine wave in the 8
% Hz LFP
th_phase = angle(hilbert(thetha));
clear MI
for i = 1:4000
Idxmodulado = find(th_phase > deg2rad(180-90) & ...
                   th_phase < deg2rad(270-90));
Idxmodulado =  Idxmodulado(randi(length(Idxmodulado),[1,i]));
% Plus some random spike indices
Idxruido = randi(length(LFP),[1,200]);
Idx = [Idxmodulado, Idxruido];

% compute MI and store it
[counts, phasebinsC] = hist(th_phase(Idx),deg2rad(-170:20:170));
p = counts/sum(counts);
MI(i) = (log(length(phasebinsC))+sum(p(p>0).*log(p(p>0))))/log(length(phasebinsC));
end

spkIn3Q = (1:4000)/4200;
subplot(111)
    plot(spkIn3Q,MI);
    xlabel('Percentage of Spikes in the Third Quadrant')
    ylabel('MI')
    axis tight
    title('MI as a Function of the Number of Spikes in the Third Quadrant')

%% Testar para modulação em função de FR in given phase bin
% manter a mesma probabilidade de um spike
% qualquer cair em um determinado fase bin. variar
% a quantidade total de spikes


%%% Toy Model Creation %%%
% setup
clear, clc, clf
f=8;
srate=1000;
dt=1/srate;
t=dt:dt:20;

% set sine waves
nonModAmp = 10;
slow_modulation_wave = sin(2*pi*7.5*t);
fast_modulated_wave = sin(2*pi*45*t).*(0.2*(slow_modulation_wave+1)+nonModAmp*0.1);

% fuse sine waves + noise to create signal
LFP = slow_modulation_wave + ...
      fast_modulated_wave  + ...
      0.24*randn(size(t));

% filter the signal for low freq for the phase
% modulating wave and for higth freq for the amp
% modulated wave
thetha = eegfilt(LFP,srate,5,10);
% gamma = eegfilt(LFP,srate,30,100);

% set index of periodic spike signal with
% consistent relation with the sine wave in the 8
% Hz LFP
th_phase = angle(hilbert(thetha));
clear MI
for j = 1:100
for i = 1:400
% Idxmodulado = find(th_phase > deg2rad(180-90) & ...
%                    th_phase < deg2rad(270-90));
% Idxmodulado =  Idxmodulado(randi(length(Idxmodulado),[1,i]));
% Plus some random spike indices
Idxruido = randi(length(LFP),[1,i]);
% Idx = [Idxmodulado, Idxruido];

% compute MI and store it
[counts, phasebinsC] = hist(th_phase(Idxruido),deg2rad(-170:20:170));
p = counts/sum(counts);
MI(j,i) = (log(length(phasebinsC))+sum(p(p>0).*log(p(p>0))))/log(length(phasebinsC));
end
end

mean_MI = mean(MI);
std_MI = std(MI);

CV = std_MI./mean_MI;

subplot(111)
    plot(1:length(MI),MI,'color',[1 1 1]/1.2), hold on
    plot(1:length(MI),mean_MI,'b-','linew',3)
    plot(1:length(MI),mean_MI+1*std_MI,'k--')
    plot(1:length(MI),mean_MI-1*std_MI,'k--'), hold off
    xlabel('Number of Spikes')
    ylabel('MI')
    axis tight
    title('Effect of the Quantity of Spikes')

% subplot(212)
%     plot(1:length(CV),CV)
%     xlabel('Number of Spikes')
%     ylabel('Coefficient of Variation')
%     axis tight
%     title('Coefficient of Variation')




% comparar com o PPC

%% Plot several von Mises distributions (varing kappa parameter)

clear, clf, clc
mu = pi;
Kappa = 0:0.2:2
for kappa = Kappa
    % compute the VonMises distribution
    phi = -pi:pi/1000:pi;
    VonMises = exp(kappa*cos(phi-mu))/(2*pi*besseli(0,kappa));
    hold on
    plot([phi+pi phi+pi+2*pi],[VonMises VonMises]*deg2rad(20))
    hold off
    set(gca,'xtick',0:pi/2:4*pi,'XTickLabel',{'0','\pi/2','\pi','3\pi/4,'})
    xlim([0 4*pi])
    xlabel('Phase (rad)')
    ylabel('Von Mises Distribution')
    title(['Von Mises distribution for varing kappa values'],'fontsize',14)
end
legend('\kappa = 0','\kappa = .2','\kappa = .4','\kappa = .6',...
    '\kappa = .8','\kappa = 1.0','\kappa = 1.2','\kappa = 1.4',...
    '\kappa = 1.6','\kappa = 1.8','\kappa = 2.0');


%% sample points from von Mises distribution (not working properly, sample biased)
% Inverse transform sampling

% given a mu and kappa values sample points from
% the respective von mises distribution and make a
% histogram
clear, clf, clc
N = 10000;
mu = pi;
kappa = 2;
phi = -pi:pi/900:pi; %d_theta precisa ser algum multiplo grande de 3
VonMises = exp(kappa*cos(phi-mu))/(2*pi*besseli(0,kappa));
CDF = cumsum(VonMises)/sum(VonMises);
phasebinsR = deg2rad(-160:20:180);
phasebinsC = deg2rad(-170:20:170);
for i=1:length(phasebinsR)
    [m, I]=min(abs(phi-phasebinsR(i)));
    phasebins_idx(i)=I;
end
% cont = 0;
% for j = 1:1000
uniformSample = sort(rand(1,N));
counts(1)=nnz(uniformSample>CDF(1) & uniformSample<=CDF(phasebins_idx(1)));
for i=2:length(phasebinsR)
        counts(i)=nnz(uniformSample>CDF(phasebins_idx(i-1)) &...
                        uniformSample<=CDF(phasebins_idx(i)));
end
amostral_probability_density = counts/sum(counts);
% aux = counts(1)>counts(18);
% cont = cont + aux;
% end
% cont
subplot(211)
    plot([phi(phasebins_idx)+pi phi(phasebins_idx)+pi+2*pi],...
        [VonMises(phasebins_idx) VonMises(phasebins_idx)]*deg2rad(20),'ko')
    set(gca,'xtick',0:pi/2:4*pi,'XTickLabel',...
        {'0','\pi/2','\pi','3\pi/4,'})
    xlim([0 4*pi])
    xlabel('Phase (rad)')
    hold off
    ylabel('value under von mises')

subplot(212)
    bar([phasebinsC+pi phasebinsC+pi+2*pi],...
        [amostral_probability_density amostral_probability_density],'k')
    set(gca,'xtick',0:pi/2:4*pi,'XTickLabel',...
        {'0','\pi/2','\pi','3\pi/4,'})
    xlim([0 4*pi])
    xlabel('Phase (rad)')
    ylabel('spk probability of each von mises value')
    
    % compute the VonMises distribution (yellow in plot)
    VonMises = exp(kappa*cos(phi-mu))/(2*pi*besseli(0,kappa));
    hold on
    plot([phi+pi phi+pi+2*pi],[VonMises VonMises]*deg2rad(20),'-y','linew',3)
    hold off
    set(gca,'xtick',0:pi/2:4*pi,'XTickLabel',{'0','\pi/2','\pi','3\pi/4,'})
    xlim([0 4*pi])
    xlabel('Phase (rad)')
    ylabel('Von Mises Distribution')
    title(['Von Mises distribution for varing kappa values'],'fontsize',14)

% dezenove = [CDF(1) CDF(phasebins_idx)]
% a = diff(dezenove)
% difA = a(9:-1:1)-a(10:18)
% VMdezenove = [VonMises(1) VonMises(phasebins_idx)]
% b = diff(VMdezenove)
% difB = b(9:-1:1)-b(10:18)
% firstH = sum(VonMises(1:900))/sum(VonMises)
% SeconH = sum(VonMises(902:end))/sum(VonMises)

%% Teste with small manual distribution to see what was the error above. Didn't found any
clear, clc
N = 10000;
chance = 0;
surrchance = 0;
for j = 1:10000
vectorTest = [6 4 2 1 2 4 6];
CDFTest = cumsum(vectorTest)/sum(vectorTest);
uniformSample = sort(rand(1,N));
counts(1)=nnz(uniformSample>0 & uniformSample<=CDFTest(1));
for i=2:length(vectorTest)
        counts(i)=nnz(uniformSample>CDFTest(i-1) &...
                        uniformSample<=CDFTest(i));
end
aux = counts(1)>counts(end);
chance = chance + aux;
aux = nnz(uniformSample<=0.5);
aux = aux>N/2;
surrchance = surrchance + aux;
end
counts
chance
surrchance
%% Draw samples from von Mises (function from online repository)

clear, clf, clc
mu = 0;
kappa = 1;
sampleSize = 20000;
phi = -pi:pi/1000:pi;
sample_points = randraw('vonmises', [mu, kappa], sampleSize);
phasebinsC = deg2rad(-170:20:170);
[counts] = hist(sample_points,phasebinsC);
amostral_probability_density = counts/sum(counts);

subplot(111)
    bar([phasebinsC+pi phasebinsC+pi+2*pi],...
        [amostral_probability_density amostral_probability_density],'k')
    set(gca,'xtick',0:pi/2:4*pi,'XTickLabel',{'0','\pi/2','\pi','3\pi/4,'})
    xlim([0 4*pi])
    xlabel('Phase (rad)')
    ylabel('spk probability')
    
    % compute the VonMises distribution (yellow in plot)
    VonMises = exp(kappa*cos(phi-mu))/(2*pi*besseli(0,kappa));
    hold on
    plot([phi+pi phi+pi+2*pi],[VonMises VonMises]*deg2rad(20),'-y','linew',3)
    hold off
    set(gca,'xtick',0:pi/2:4*pi,'XTickLabel',{'0','\pi/2','\pi','3\pi/4,'})
    xlim([0 4*pi])
    xlabel('Phase (rad)')
    ylabel('Von Mises Distribution')
    title(['Von Mises distribution for varing kappa values'],'fontsize',14)


%% PLV for diferent numbers of observations

clear,clf,clc
Kappa = [0.1 1/2 1 2 4 8];
count = 0;
mu = 0;
sampleSize = 100;
color1 = [0*(1:6)' (.9-.1*(1:6))' 0*(1:6)'];
color2 = [0*(1:6)' 0*(1:6)' (.9-.1*(1:6))'];
color3 = [(.9-.1*(1:6))' 0*(1:6)' 0*(1:6)'];
tic
clear PLV PPCAll
for kappa = Kappa
count = count + 1
for trial = 1:10
    sample_points = randraw('vonmises', [mu, kappa], sampleSize) + pi;
    for n = 1:sampleSize
        % PLV %
        PLV(count,n,trial) = abs(mean(exp(1i*sample_points(1:n))));
        % PPC %
        PPC = 0;
        for j = 1:n
            for k = j+1:n
                if j~=k
                    PPC = PPC + real(exp(1i*(sample_points(j)-sample_points(k))));
                end
            end
        end
        PPC = PPC*2/(n*(n-1));
        PPCAll(count,n,trial) = PPC;
        % Converting phase values to the polar representation
        crss = exp(1i*(sample_points));
        % Getting the number of vectors
        dof = sum(~isnan(crss));
        % Magnitude of the imaginary part
        sinSum = abs(nansum(imag(crss)));
        % Magnitude of the real part
        cosSum = nansum(real(crss));
        % PPC calculation
        ppc = (cosSum.^2+sinSum.^2 - dof)./(dof.*(dof-1));
        ppcAll(count,n,trial) = ppc;
    end
end
end
meanPLV = mean(PLV,3);
meanPPCAll = mean(PPCAll,3);
meanppcAll = mean(ppcAll,3);
toc

subplot(111), clf
    for i=1:length(Kappa)
        hold on
        plot(1:sampleSize,meanPLV(i,:),'color',color1(i,:),'linew',2)
        plot(1:sampleSize,meanPPCAll(i,:),'color',color2(i,:),'linew',2)
        plot(1:sampleSize,meanppcAll(i,:),'color',color3(i,:),'linew',2)
    end
    hold off

    % comparar com dados reais e verificar se a
    % quantidade de 
    
    % procurar comparações do PPC
    
    
%%
% setup
clear, clc, clf
load('.\Data.mat\SpkBuz.mat')
load('.\Data.mat\thetaBuz.mat')
load('.\Data.mat\LFPBuz.mat')
LFP = LFP';

for lo = 1:40
lo
% Filter LFP to wanted frequency
LFP_filt = eegfilt(LFP,srate,lo,lo+4);

% Extract phases
phases = angle(hilbert(LFP_filt));

% get the spike times and there indexes
spktimes = Raster{7}; % in s
spkind = round(spktimes*srate);  % round ms to indexes

% count the number of spikes in each bin
phasebins = deg2rad(-170:20:170);
[counts] = hist(phases(spkind),phasebins);

% scale dow the count to get a prob. distribution
p = counts/sum(counts);

% get kappa with tool kit 
[mu,PLV(lo),sigma,CI,kappa(lo)]=anglemean(phases(spkind));

% compute MI and store it
MI(lo) = (log(length(phasebins))+sum(p(p>0).*log(p(p>0)))) ...
     /log(length(phasebins));

% % compute the VonMises distribution (yellow in plot)
% phi = -pi:pi/1000:pi;
% VonMises = exp(kappa*cos(phi-mu))/(2*pi*besseli(0,kappa));

% Converting phase values to the polar representation
crss = exp(1i*(phases));
% Getting the number of vectors
dof = sum(~isnan(crss));
% Magnitude of the imaginary part
sinSum = abs(nansum(imag(crss)));
% Magnitude of the real part
cosSum = nansum(real(crss));
% PPC calculation
ppc(lo) = (cosSum.^2+sinSum.^2 - dof)./(dof.*(dof-1));
end

subplot(411)
    plot((1:40)+2,PLV)
    xlabel('Freq (Hz)')
    ylabel('PLV','fontsize',20)
    title(['PLV_{max} = ' num2str(max(PLV))])

subplot(412)
    plot((1:40)+2,kappa)
    xlabel('Freq (Hz)')
    ylabel('\kappa','fontsize',20)
    title(['\kappa_{max} = ' num2str(max(kappa))])

subplot(413)
    plot((1:40)+2,MI)
    xlabel('Freq (Hz)')
    ylabel('MI','fontsize',20)
    title(['MI_{max} = ' num2str(max(MI))])

% One can see that MI goes linear with PLV^2 and kappa^2 
subplot(414)
    plot((1:40),ppc);
    xlabel('Freq (Hz)')
    ylabel('PPC','fontsize',20)
    title(['MI_{max} = ' num2str(max(ppc))])

subplot(111)
    % bar(phasebins,counts,'k')
    bar([phasebins+pi phasebins+pi+2*pi],...
        [p p],'k')
    hold on
    plot([phi+pi phi+pi+2*pi],[VonMises VonMises]*deg2rad(20),'y-','linew',3)
    hold off
    set(gca,'xtick',0:pi/2:4*pi,'XTickLabel',{'0','\pi/2','\pi','3\pi/4,'})
    xlim([0 4*pi])
    xlabel('Phase (rad)')
    ylabel('spk probability')
    title(['Kappa_{Von Mises} = ' num2str(kappa)],'fontsize',14)
























