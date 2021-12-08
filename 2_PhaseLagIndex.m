%% Segundo Mes Stam
% Use surrogate data to compute a z-score
% Find out how to avoid DeltaPhase jumps on randn data
% Think on testes to do in the data, what do i whant to show?
% Onde entra o efeito do sample size? Se aumentar a janela
% menor a chance de conseguir uma diferença de fase bounded


%% PLI Works

% Setup
clear, clc
srate = 1000; f = 8; % in Hz
dt = 1/srate; Tmax = 01; t = dt:dt:Tmax; % in Seconds
N = length(t);

% Signal
phaseLag = (pi*randn);
X = sin(2*pi*f.*t);
Y = sin(2*pi*f.*t + phaseLag);

% Analytical Signal (ARepresentation)
RAx = hilbert(X);
RAy = hilbert(Y);

% Hilbert Transform part
HTx = imag(RAx);
HTy = imag(RAy);

% Real time Angle
Angx = angle(RAx);
Angy = angle(RAy);

% Real time amplitude
Ampx = abs(RAx);
Ampy = abs(RAy);

% Phase Syncronization
% DeltaPhase = unwrap(Angx)-unwrap(Angy); % não recomendado
DeltaPhase = angle(exp(1i*(Angx-Angy))); % recomendado

% R index
R = mean(DeltaPhase);
disp(rad2deg(abs(R)));
disp(rad2deg(abs(phaseLag)));

% plot([real(exp(1i*R))],[imag(exp(1i*R))],'ro','markerf','r')
% hold on
% plot([0 real(exp(1i*R))],[0 imag(exp(1i*R))],'r-')
% plot([real(exp(1i*phaseLag))],[imag(exp(1i*phaseLag))],'bo','markerf','b')
% plot([0 real(exp(1i*phaseLag))],[0 imag(exp(1i*phaseLag))],'b-') 
% plot(exp(1i*(0:0.01:2*pi)),'k-')
% axis square

%%% PLI %usar a função sign()
if(DeltaPhase>-pi & DeltaPhase<pi)
    disp('PLI calculated with -pi<DeltaPhase<pi')
    above0 = sum(DeltaPhase>0);
    bellow0 = N-above0;
    PLI = abs((above0-bellow0)/N);
elseif (DeltaPhase>0 & DeltaPhase<2*pi)
    disp('PLI calculated with 0<DeltaPhase<2*pi')
    above0 = sum(sin(DeltaPhase)>0);
    bellow0 = N-above0;
    PLI = abs((above0-bellow0)/N);
else
    PLI = 0;
    disp('PLI error: DeltaPhase <-pi ou >2*pi')
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % computing the Phase-Locking Value
    PLV = abs(mean(exp(1i*DeltaPhase)));
    disp(['PLV = ' num2str(PLV)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=figure(1); clf; set(h,'WindowStyle','docked')
subplot(2,1,1)
    plot(t,real(RAx))
    hold on
    plot(t,real(RAy))
    plot(t,HTx,'g-')
    plot(t,HTy,'r-')
    plot(t,Ampx,'k-','linewidth',2)
    plot(t,Ampy,'k-','linewidth',2)
    xlabel('Time (s)')
    ylabel('mV')
    title(['Signals with Phase Lag of ' num2str(phaseLag) ' rad'])
    hold off
subplot(2,1,2)
    plot(t,Angx)
    hold on
    plot(t,Angy)
    plot(t,DeltaPhase,'yo-')
    xlabel('Time (s)')
    ylim([-2*pi 2*pi])
    ylabel('\Phi (rad)')
    set(gca,'ytick',-pi:pi/2:pi,...
        'yticklabel',{'-pi','-pi/2','0','pi/2','pi'})
    title(['PLI = ' num2str(mean(PLI))])
    hold off
    
% sync = co_sync(Angx,Angy,1,1);

%% (PLI janelado) >>> não usamos janelamento no calculo de PLI
% janelamento não é usado em calculo de PLI
% PLI per win
win = 1*srate; % window for PLI calculation
step = 0.1*win;
Nwin = (length(t)-win)/step+1;
T = zeros(1,Nwin);
for nwin = 1:Nwin
nwin;
    winidx = (1:win) + (nwin-1)*step; 
    N = length(winidx);
    if(DeltaPhase(winidx)>-pi & DeltaPhase(winidx)<pi)
        disp('east')
        above0 = sum(DeltaPhase(winidx)>0);
        bellow0 = N-above0;
        PLI(nwin) = abs((above0-bellow0)/N);
    elseif (DeltaPhase(winidx)>0 & DeltaPhase(winidx)<2*pi)
        disp('west')
        above0 = sum(sin(DeltaPhase(winidx))>0);
        bellow0 = N-above0;
        PLI(nwin) = abs((above0-bellow0)/N);
    else
        PLI(nwin) = 0;
        disp('PLI error')
    end
    T(nwin) = nwin; 
end

subplot(2,1,2)
    plot(T(3:19),PLI(3:19))
    ylim([0 1])
%% Filtered white noise 
% gives low PLI and PLV as expected

% Setup
clear, clc
srate = 1000; % in Hz
dt = 1/srate; Tmax = 3; t = dt:dt:Tmax; % in Seconds
N = length(t);

% Signal
X = randn(size(t));
X = eegfilt(X,srate,7,9);
Y = randn(size(t));
Y = eegfilt(Y,srate,7,9);

RAx = hilbert(X);   RAy = hilbert(Y);
HTx = imag(RAx);    HTy = imag(RAy);
Angx = angle(RAx);  Angy = angle(RAy);
Ampx = abs(RAx);    Ampy = abs(RAy);

DeltaPhase = angle(exp(1i*(Angx-Angy))); % recomendado
R = mean(DeltaPhase);
disp(rad2deg(abs(R)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% computing the Phase-Lag Index %%%
if(DeltaPhase>-pi & DeltaPhase<pi)
    disp('PLI calculated with -pi<DeltaPhase<pi')
    above0 = sum(DeltaPhase>0);
    bellow0 = N-above0;
    PLI = abs((above0-bellow0)/N);
elseif (DeltaPhase>0 & DeltaPhase<2*pi)
    disp('PLI calculated with 0<DeltaPhase<2*pi')
    above0 = sum(sin(DeltaPhase)>0);
    bellow0 = N-above0;
    PLI = abs((above0-bellow0)/N);
else
    PLI = 0;
    disp('PLI error: DeltaPhase <-pi ou >2*pi')
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % computing the Phase-Locking Value
    PLV = abs(mean(exp(1i*DeltaPhase)));
    disp(['PLV = ' num2str(PLV)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=figure(1); clf; set(h,'WindowStyle','docked')
subplot(2,1,1)
    plot(t,real(RAx))
    hold on
    plot(t,real(RAy))
    plot(t,HTx,'g-')
    plot(t,HTy,'r-')
    plot(t,Ampx,'k-','linewidth',2)
    plot(t,Ampy,'k-','linewidth',2)
    xlabel('Time (s)')
    ylabel('mV')
    title('Filtered White Noise')
    hold off
subplot(2,1,2)
    plot(t,Angx)
    hold on
    plot(t,Angy)
    xlabel('Time (s)')
    ylim([-2*pi 2*pi])
    ylabel('\Phi (rad)')
    set(gca,'ytick',-pi:pi/2:pi,...
        'yticklabel',{'-pi','-pi/2','0','pi/2','pi'})
    title(['PLI = ' num2str(PLI)])
    plot(t,DeltaPhase,'yo-')
    hold off

%% Windowled Real Data >>> Pode ser afetado pelo bug da fase?
%verificar qual tamanho de janela é melhor para cada situação


% Setup
clear, clc, load('.\Data.mat\LFPprobe.mat')
dt = 1/srate; Tmax = 4; t = dt:dt:Tmax; % in Seconds
N = length(t);
ff = 9; nyquist = srate/2;

% Signal
X = eegfilt(double(LFPprobe(1,1:N)),srate,8,10);
Y = eegfilt(double(LFPprobe(16,1:N)),srate,8,10);

% Kernell = sin(2*pi*ff*t(1:500));
% X = conv(LFPprobe(1, 1:N), Kernell, 'same');
% Y = conv(LFPprobe(16,1:N), Kernell, 'same');

RAx = hilbert(X);   RAy = hilbert(Y);
Angx = angle(RAx);  Angy = angle(RAy);

DeltaPhase = angle(exp(1i*(Angx-Angy))); % recomendado
R = mean(DeltaPhase); disp(rad2deg(abs(R)));

% PLI per win
win = 1*srate; % window for PLI calculation
step = 0.2*win;
Nwin = (length(t)-win)/step+1;
T = zeros(1,Nwin); 
for nwin = 1:Nwin
nwin;
    winidx = (1:win) + (nwin-1)*step;
    N = length(winidx);
    if(DeltaPhase(winidx)>-pi & DeltaPhase(winidx)<pi)
        disp('PLI calculated with -pi<DeltaPhase<pi')
        above0 = sum(DeltaPhase(winidx)>0);
        bellow0 = N-above0;
        PLI(nwin) = abs((above0-bellow0)/N);
    elseif (DeltaPhase(winidx)>0 & DeltaPhase(winidx)<2*pi)
        disp('PLI calculated with 0<DeltaPhase<2*pi')
        above0 = sum(sin(DeltaPhase(winidx))>0);
        bellow0 = N-above0;
        PLI(nwin) = abs((above0-bellow0)/N);
    else
        PLI(nwin) = 0;
        disp('PLI error: DeltaPhase <-pi ou >2*pi')
    end
    T(nwin) = nwin; 
end

h=figure(1); clf; set(h,'WindowStyle','docked')
subplot(3,1,1)
    plot(t,real(RAx))
    hold on
    plot(t,real(RAy))
    xlabel('Time (s)')
    ylabel('mV')
    title('9 Hz Filtered Real Data')
    hold off
    
subplot(3,1,2)    
    plot(t,Angx)
    hold on
    plot(t,Angy)
    xlabel('Time (s)')
    ylim([-2*pi 2*pi])
    ylabel('\Phi (rad)')
    set(gca,'ytick',-pi:pi/2:pi,...
        'yticklabel',{'-pi','-pi/2','0','pi/2','pi'})
    title(['Mean PLI = ' num2str(mean(PLI))])
    plot(t,DeltaPhase,'yo-')
    hold off
    
subplot(3,1,3)
    plot(T,PLI)
    xlabel('Window')
    ylabel('PLI')
    ylim([0 1])
    title('Windowled Real Data')

%% noisy sine waves >>> Sensible to Noise

% Setup
clear, clc
srate = 1000; f = 8; % in Hz
dt = 1/srate; Tmax = 1; t = dt:dt:Tmax; % in Seconds
N = length(t);

% Signal
phaseLag = pi/2; % MUDAR O ANGULO ESTÁ MUDANDO O PLI
X = sin(2*pi*f.*t)+0.3*randn(size(t));
Y = sin(2*pi*f.*t + phaseLag+0.3*randn(size(t)));

RAx = hilbert(X);   RAy = hilbert(Y);
HTx = imag(RAx);    HTy = imag(RAy);
Angx = angle(RAx);  Angy = angle(RAy);

DeltaPhase = angle(exp(1i*(Angx-Angy))); % recomendado

% R index
R = mean(DeltaPhase);
disp(rad2deg(abs(R)));
disp(rad2deg(abs(phaseLag)));

%%% PLI
if(DeltaPhase>-pi & DeltaPhase<pi)
    disp('PLI calculated with -pi<DeltaPhase<pi')
    above0 = sum(DeltaPhase>0);
    bellow0 = N-above0;
    PLI = abs((above0-bellow0)/N);
elseif (DeltaPhase>0 & DeltaPhase<2*pi)
    disp('PLI calculated with 0<DeltaPhase<2*pi')
    above0 = sum(sin(DeltaPhase)>0);
    bellow0 = N-above0;
    PLI = abs((above0-bellow0)/N);
else
    PLI = 0;
    disp('PLI error: DeltaPhase <-pi ou >2*pi')
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % computing the Phase-Locking Value
    PLV = abs(mean(exp(1i*DeltaPhase)));
    disp(['PLV = ' num2str(PLV)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=figure(1); clf; set(h,'WindowStyle','docked')
subplot(2,1,1)
    plot(t,real(RAx))
    hold on
    plot(t,real(RAy))
    plot(t,HTx,'g-')
    plot(t,HTy,'r-')
    xlabel('Time (s)')
    ylabel('mV')
    title('Noisy Sine Waves')
    hold off
subplot(2,1,2)
    plot(t,Angx)
    hold on
    plot(t,Angy)
    xlabel('Time (s)')
    ylim([-2*pi 2*pi])
    ylabel('\Phi (rad)')
    set(gca,'ytick',-pi:pi/2:pi,...
        'yticklabel',{'-pi','-pi/2','0','pi/2','pi'})
    title(['PLI = ' num2str(PLI)])
    plot(t,DeltaPhase,'yo-')
    hold off

%% Finite sample bias >> visualization

clear, clc, clf
srate = 1000; f = 8;
dt = 1/srate; t = dt:dt:2;

ruido1 = 1*randn(size(t));
ruido2 = 1*randn(size(t));

LFP1 = sin(2*pi*f*t) + ruido1;
LFP2 = sin(2*pi*f*t+pi/2) + ruido2;

LFP1filtrado = eegfilt(LFP1,srate,60,70);
LFP2filtrado = eegfilt(LFP2,srate,60,70);

% LengthVector = 0.005:0.005:4;
LengthVector = 5:5:2000;
count = 0;
PLV_epochlength = size(LengthVector);

for EpochLength = LengthVector; % in seconds
    count = count+1;

    Phase1 = angle(hilbert(LFP1filtrado(1:EpochLength)));
    Phase2 = angle(hilbert(LFP2filtrado(1:EpochLength)));

    DeltaPhase = angle(exp(1i*(Phase2-Phase1))); % metodo recomendado
    PLV = abs(mean(exp(1i*DeltaPhase)));
    PLV_epochlength(count) = PLV;
end

plot(LengthVector,PLV_epochlength)

%% Finite sample bias >> quntification by variance

clear, clc, clf
srate = 1000; f = 8;
dt = 1/srate; t = dt:dt:2;

N = 100;
% LengthVector = 0.005:0.005:4;
LengthVector = 5:5:2000;
PLVspec_holder = zeros(N, length(LengthVector));

for i = 1:N
ruido1 = 1*randn(size(t));
ruido2 = 1*randn(size(t));

LFP1 = sin(2*pi*f*t) + ruido1;
LFP2 = sin(2*pi*f*t+pi/2) + ruido2;

LFP1filtrado = eegfilt(LFP1,srate,60,70);
LFP2filtrado = eegfilt(LFP2,srate,60,70);

count = 0;
PLV_epochlength = size(LengthVector);

for EpochLength = LengthVector; % in seconds
    count = count+1;

    Phase1 = angle(hilbert(LFP1filtrado(1:EpochLength)));
    Phase2 = angle(hilbert(LFP2filtrado(1:EpochLength)));

    DeltaPhase = angle(exp(1i*(Phase2-Phase1))); % metodo recomendado
    PLV = abs(mean(exp(1i*DeltaPhase)));
    PLV_epochlength(count) = PLV;
end
% plot(LengthVector,PLV_epochlength)
PLVspec_holder(i,:) = PLV_epochlength;
end
V = var(PLVspec_holder);

plot(1:length(LengthVector),V);
























