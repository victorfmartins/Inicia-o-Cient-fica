%% Segundo Mes Vink 

clear, clc
srate = 1000;
dt = 1/srate; Tmax = 1; t = dt:dt:Tmax; % in Seconds

Ntrials = 10; % j is a trial
Tsamples = 1*srate; % per trial
Sensors = 2; % serie temporal de 2 sensores simultaneos por trial
Ksources = 10; % fontes a serem misturadas para formar os sinais.

% S has Ntrials pairs of signals (Sensors) with Tsamples points each
% Sj is the matrix S(1,:,:) of size 2xT
S = zeros(Ntrials,Sensors,Tsamples); 
V = zeros(Ntrials,Ksources,Tsamples);
% Vj = zeros(Ksources,Tsamples); % matrix that model each signal

% modificar as sources para ter trials diferentes
for j = 1:Ntrials
    for k = 1:Ksources
        V(j,k,:) = sin(2*pi*(k+10*j-10)*t)+0.1*randn(1,Tsamples);
    end
end

figure(1), clf
    j = 1; % choose trial to plot
    plot(t,squeeze(V(j,1,:)), 'linew', 1); hold on
    plot(t,squeeze(V(j,2,:)), 'linew', 2)
    plot(t,squeeze(V(j,3,:)), 'linew', 3)
    plot(t,squeeze(V(j,4,:)), 'linew', 4)
    legend('1','2','3','4')
    xlabel('Time (s)')
    ylabel('Power (\muV)')
    title('Visualization of First 4 Sources of First Trial')

A = zeros(Sensors,Ksources);
% A = randn(Sensors,Ksources);
for k = 1:Ksources
    if mod(k,2)
        A(1,k) = 1;
    else
        A(2,k) = 1;
    end
end

for j = 1:Ntrials
    S(j,:,:) = A*squeeze(V(j,:,:));
end

figure(2), clf
    subplot(2,1,1)
    j = 1; % choose trial to plot
    a = [[1     0     1     0     1     0     1     0     1     0];...
         [0     1     0     1     0     1     0     1     0     1]];
    Sj = a*squeeze(V(j,:,:));
    plot(t, Sj(1,:))
    subplot(2,1,2)
    Sj4visualisation = sin(2*pi*(1+j-1)*t)+sin(2*pi*(3+j-1)*t)...
        +sin(2*pi*(5+j-1)*t)+sin(2*pi*(7+j-1)*t)+sin(2*pi*(9+j-1)*t);
    plot(t,Sj4visualisation)


Y = zeros(Ntrials,Ksources,Tsamples); % fft(sources)
% Yj ? (Y(1,j),…,Y(K,j))' is a vector of 
% complex-valued Fourier spectra whose 
% k-th element is obtained by DFT-ing 
% the k-th sensor row of Vj
% Yj = zeros(Ksources,Tsamples); %Yj(1:K,j)' = cada pt é um F. Spectr
% fft trabalha na coluna. A coluna precisa 
% ser uma time series do mesmo sinal, então 
% vamos tomar a transposta de Vj.
% A fft(Vj') retorna 1000 linhas de freq,
% das quais apenas 500 são importantes, e
% 10 colunas.
for j = 1:Ntrials
% modificar as sources para ter trials diferentes
% Yj(:,:) = fft(Vj')'; 
Y(j,:,:) = fft(squeeze(V(j,:,:))')'/length(t);
end

figure(3), clf
    F = linspace(0,srate/2,floor(length(t)/2)+1); 
    subplot(4,1,1)
    plot(F,abs(2*squeeze(Y(1,1,1:length(F)))), 'linew', 1); hold on
    y=fft(sin(2*pi*1*t))/length(t); 
    plot(F,2*abs(y(1:length(F))), 'r', 'linew', 1); hold off
    xlim([0 10])
    xlabel('Frequency (Hz)'), ylabel('Power (\muV)')
    title('Frequency of First 4 Sources of First Trial')
    subplot(4,1,2)
    plot(F,abs(2*squeeze(Y(1,2,1:length(F)))), 'linew', 1); hold on
    y=fft(sin(2*pi*2*t))/length(t); 
    plot(F,2*abs(y(1:length(F))), 'r', 'linew', 1); hold off
    xlim([0 10])
    xlabel('Frequency (Hz)'), ylabel('Power (\muV)')
    subplot(4,1,3)
    plot(F,abs(2*squeeze(Y(1,3,1:length(F)))), 'linew', 1); hold on
    y=fft(sin(2*pi*3*t))/length(t); 
    plot(F,2*abs(y(1:length(F))), 'r', 'linew', 1); hold off
    xlim([0 10])
    xlabel('Frequency (Hz)'), ylabel('Power (\muV)')
    subplot(4,1,4)
    plot(F,abs(2*squeeze(Y(1,4,1:length(F)))), 'linew', 1); hold on
    y=fft(sin(2*pi*4*t))/length(t); 
    plot(F,2*abs(y(1:length(F))), 'r', 'linew', 1); hold off
    xlim([0 10]), ylim([0 1])
    xlabel('Frequency (Hz)'), ylabel('Power (\muV)')

% Zj = zeros(Sensors,Tsamples); % prealocation not recommended
Z = zeros(Ntrials,Sensors,Tsamples);

for j = 1:Ntrials
Zj = A*squeeze(Y(j,:,:)); % Zj = fft(Sj)
Z(j,:,:) = Zj; % N x 2 x T
end

figure(4), clf % FORMA ALTERNATIVA: Zj = fft(Sj)
    subplot(2,1,1)
    Zj = A*squeeze(Y(1,:,:));
    plot(F, 2*abs(squeeze(Z(1,1,1:length(F))))); hold on
    plot(F, 2*abs(squeeze(Z(1,2,1:length(F)))),'r'); hold off
    % plot(F,2*abs(Zj(1,1:length(F))),'ko-','linewidth',3,...
    %     'markersize',14,'markerfacecolor','w')
    % stem(F,2*abs(Zj(1,1:length(F))),'k-','linew',3)
    xlim([0 20]), ylim([0 1])
    xlabel('Frequency (Hz)'), ylabel('Power (\muV)')
    title('Spectrum via Source')
    subplot(2,1,2)
    Zj4visualisation = fft(squeeze(S(1,1,:)))/length(t);
    plot(F,2*abs(Zj4visualisation(1:length(F)))); hold on
    Zj4visualisation = fft(squeeze(S(1,2,:)))/length(t);
    plot(F,2*abs(Zj4visualisation(1:length(F))),'r'); hold off
    xlim([0 20]), ylim([0 1])
    xlabel('Frequency (Hz)'), ylabel('Power (\muV)')
    title('Spectrum via Signal')

X = zeros(Ntrials, 2*Tsamples-1);
% non-diagonal part of the cross-spectrum                   ???????
for j = 1:Ntrials
Zconj = conj(squeeze(Z(j,2,:)));
X(j,:) = xcorr(squeeze(Z(j,1,:)),Zconj); % X=R*exp(i*theta),
end

figure(5), clf
    subplot(2,1,1)
    plot(F, 2*abs(squeeze(Z(1,1,1:length(F))))); hold on
    plot(F, 2*abs(squeeze(Z(1,2,1:length(F))))); hold off
    xlim([0 20])
    subplot(2,1,2)
    zconj = conj(squeeze(Z(1,2,:)));
    [x, lags] = xcorr(squeeze(Z(1,1,:)),zconj);
    plot(lags,x,'ko-')
    % plot3(lags,real(x), imag(x),'ko-');
    % polar(angle(x), abs(x))
    xlim([0 20])

% If two variables are uncorrelated, 
% there is no linear relationship between them
% For (linearly) uncorrelated source activities, we have E{YkYl?}=0
Y1 = squeeze(Y(1,:,:));
Y2 = squeeze(Y(2,:,:));
uncorr = Y1*Y2';
meancorr = mean(uncorr); %should yeild zero
soma = sum(uncorr);

figure(6), clf
    plot(1:10,real(meancorr)), hold on
    plot(1:10,real(soma)), hold off
    legend('mean', 'sum')

% Coherency
C = zeros(Ntrials,1);
for j=1:Ntrials
Cj = mean(X(j,:))/(mean(abs(squeeze(Z(j,1,:)).^2))*...
              mean(abs(squeeze(Z(j,2,:)).^2)));
C(j) = Cj;
end

% Phase Lock Value (PLV)
PLV = zeros(Ntrials,1);
for j=1:Ntrials
PLVj = mean(exp(1i*angle(X(j,:)))); % usar diff(angle)?
PLV(j) = PLVj;
end

% PLI
PLI = zeros(Ntrials,1);
for j=1:Ntrials
PLIj = abs(mean(sign(imag(X(j,:)))));
PLI(j) = PLIj;
end

% WPLI
WPLI = zeros(Ntrials,1);
for j=1:Ntrials
WPLIj = abs(mean(imag(X(j,:))))/mean(abs(imag(X(j,:)))); 
WPLI(j) = WPLIj;
end
% WPLI = abs(mean(imsg(X))*

%% Sumarization:

clear, clc
srate = 1000;
dt = 1/srate; Tmax = 1; t = dt:dt:Tmax; % in Seconds

Ntrials = 10; % j is a trial
Tsamples = 1*srate; % per trial
Sensors = 2; % serie temporal de 2 sensores simultaneos por trial
Ksources = 2; % fontes a serem misturadas para formar os sinais.

% Model Declarations
A = [[1 0];[0 1]]; % Sensors x Ksources
S = zeros(Ntrials,Sensors,Tsamples); % Ntrials pairs = signal;
V = zeros(Ntrials,Ksources,Tsamples); % Ntrials Ksources each = model;
Y = zeros(Ntrials,Ksources,Tsamples); % Y = fft(V)
Z = zeros(Ntrials,Sensors,Tsamples); % Z = fft(S)
X = zeros(Ntrials, 2*Tsamples-1); % X = xcorr(Z)

% Calculations
for j = 1:Ntrials
    for k = 1:Ksources
        V(j,k,:) = sin(2*pi*(k+10*j-10)*t)+0.1*randn(1,Tsamples);
    end
    S(j,:,:) = A*squeeze(V(j,:,:));
    Y(j,:,:) = fft(squeeze(V(j,:,:))')'/length(t);
    Zj = A*squeeze(Y(j,:,:)); % Zj = fft(Sj)
    Z(j,:,:) = Zj; % N x 2 x T
    Zconj = conj(squeeze(Z(j,2,:)));
    X(j,:) = xcorr(squeeze(Z(j,1,:)),Zconj);
end

% Metrics Declarations
C = zeros(Ntrials,1);
PLV = zeros(Ntrials,1);
PLI = zeros(Ntrials,1);
WPLI = zeros(Ntrials,1);

% Calculations
for j=1:Ntrials
    Cj = mean(X(j,:))/(mean(abs(squeeze(Z(j,1,:)).^2))*...
                       mean(abs(squeeze(Z(j,2,:)).^2)));
    PLVj = mean(exp(1i*angle(X(j,:)))); % usar diff(angle)?
    PLIj = abs(mean(sign(imag(X(j,:)))));
    WPLIj = abs(mean(imag(X(j,:))))/mean(abs(imag(X(j,:)))); 
    
    C(j) = Cj;
    PLV(j) = PLVj;
    PLI(j) = PLIj;
    WPLI(j) = WPLIj;
end

Y1 = squeeze(Y(1,:,:));
Y2 = squeeze(Y(2,:,:));
uncorr = Y1*Y2';
meancorr = mean(uncorr); %should yeild zero
soma = sum(uncorr);
%{
figure(1), clf
    j = 1; % choose trial to plot
    plot(t,squeeze(V(j,1,:)), 'linew', 1); hold on
    plot(t,squeeze(V(j,2,:)), 'linew', 2)
    legend('1','2')
    xlabel('Time (s)')
    ylabel('Power (\muV)')
    title('Visualization of First 4 Sources of First Trial')

figure(2), clf
    subplot(2,1,1)
    j = 1; % choose trial to plot
    a = [[1     0];...
         [0     1]];
    Sj = a*squeeze(V(j,:,:));
    plot(t, Sj(1,:))
    subplot(2,1,2)
    Sj4visualisation = sin(2*pi*(1+j-1)*t);
    plot(t,Sj4visualisation)

figure(3), clf
    F = linspace(0,srate/2,floor(length(t)/2)+1); 
    subplot(2,1,1)
    plot(F,abs(2*squeeze(Y(1,1,1:length(F)))), 'linew', 1); hold on
    y=fft(sin(2*pi*1*t))/length(t); 
    plot(F,2*abs(y(1:length(F))), 'r', 'linew', 1); hold off
    xlim([0 10])
    xlabel('Frequency (Hz)'), ylabel('Power (\muV)')
    title('Frequency of First 4 Sources of First Trial')
    subplot(2,1,2)
    plot(F,abs(2*squeeze(Y(1,2,1:length(F)))), 'linew', 1); hold on
    y=fft(sin(2*pi*2*t))/length(t); 
    plot(F,2*abs(y(1:length(F))), 'r', 'linew', 1); hold off
    xlim([0 10])
    xlabel('Frequency (Hz)'), ylabel('Power (\muV)')
    
figure(4), clf % FORMA ALTERNATIVA: Zj = fft(Sj)
    subplot(2,1,1)
    Zj = A*squeeze(Y(1,:,:));
    plot(F, 2*abs(squeeze(Z(1,1,1:length(F))))); hold on
    plot(F, 2*abs(squeeze(Z(1,2,1:length(F)))),'r'); hold off
    % plot(F,2*abs(Zj(1,1:length(F))),'ko-','linewidth',3,...
    %     'markersize',14,'markerfacecolor','w')
    % stem(F,2*abs(Zj(1,1:length(F))),'k-','linew',3)
    xlim([0 20]), ylim([0 1])
    xlabel('Frequency (Hz)'), ylabel('Power (\muV)')
    title('Spectrum via Source')
    subplot(2,1,2)
    Zj4visualisation = fft(squeeze(S(1,1,:)))/length(t);
    plot(F,2*abs(Zj4visualisation(1:length(F)))); hold on
    Zj4visualisation = fft(squeeze(S(1,2,:)))/length(t);
    plot(F,2*abs(Zj4visualisation(1:length(F))),'r'); hold off
    xlim([0 20]), ylim([0 1])
    xlabel('Frequency (Hz)'), ylabel('Power (\muV)')
    title('Spectrum via Signal')

figure(5), clf
    subplot(2,1,1)
    plot(F, 2*abs(squeeze(Z(1,1,1:length(F))))); hold on
    plot(F, 2*abs(squeeze(Z(1,2,1:length(F))))); hold off
    xlim([0 20])
    subplot(2,1,2)
    zconj = conj(squeeze(Z(1,2,:)));
    [x, lags] = xcorr(squeeze(Z(1,1,:)),zconj);
    plot(lags,x,'ko-')
    % plot3(lags,real(x), imag(x),'ko-');
    % polar(angle(x), abs(x))
    xlim([0 20])

figure(6), clf
    plot(1:2,real(meancorr)), hold on
    plot(1:2,real(soma)), hold off
    legend('mean', 'sum')
%}

%% Variando Ntrials

clear, clc
srate = 1000;
dt = 1/srate; Tmax = 1; t = dt:dt:Tmax; % in Seconds

Ntrials = 10; % j is a trial
Tsamples = 1*srate; % per trial
Sensors = 2; % serie temporal de 2 sensores simultaneos por trial
Ksources = 2; % fontes a serem misturadas para formar os sinais.

% Model Declarations
A = [[1 0];[0 1]]; % Sensors x Ksources
S = zeros(Ntrials,Sensors,Tsamples); % Ntrials pairs = signal;
V = zeros(Ntrials,Ksources,Tsamples); % Ntrials Ksources each = model;
Y = zeros(Ntrials,Ksources,Tsamples); % Y = fft(V)
Z = zeros(Ntrials,Sensors,Tsamples); % Z = fft(S)
X = zeros(Ntrials, 2*Tsamples-1); % X = xcorr(Z)

% Calculations
for ntrials = 2:Ntrials
for j = 1:ntrials
    for k = 1:Ksources
        V(j,k,:) = sin(2*pi*(k+10*j-10)*t)+0.1*randn(1,Tsamples);
    end
    S(j,:,:) = A*squeeze(V(j,:,:));
    Y(j,:,:) = fft(squeeze(V(j,:,:))')'/length(t);
    Zj = A*squeeze(Y(j,:,:)); % Zj = fft(Sj)
    Z(j,:,:) = Zj; % N x 2 x T
    Zconj = conj(squeeze(Z(j,2,:)));
    X(j,:) = xcorr(squeeze(Z(j,1,:)),Zconj);
end

% Metrics Declarations
C = zeros(Ntrials,1);
PLV = zeros(Ntrials,1);
PLI = zeros(Ntrials,1);
WPLI = zeros(Ntrials,1);

% Calculations
for j=1:Ntrials
    Cj = mean(X(j,:))/(mean(abs(squeeze(Z(j,1,:)).^2))*...
                       mean(abs(squeeze(Z(j,2,:)).^2)));
    PLVj = mean(exp(1i*angle(X(j,:)))); % usar diff(angle)?
    PLIj = abs(mean(sign(imag(X(j,:)))));
    WPLIj = abs(mean(imag(X(j,:))))/mean(abs(imag(X(j,:)))); 
    
    C(j) = Cj;
    PLV(j) = PLVj;
    PLI(j) = PLIj;
    WPLI(j) = WPLIj;
end

% Plots
%plot

end

%% Variando Tsamples

clear, clc
srate = 1000;
dt = 1/srate; Tmax = 1; t = dt:dt:Tmax; % in Seconds

Ntrials = 10; % j is a trial
Tsamples = 1*srate; % per trial
Sensors = 2; % serie temporal de 2 sensores simultaneos por trial
Ksources = 2; % fontes a serem misturadas para formar os sinais.

% Model Declarations
A = [[1 0];[0 1]]; % Sensors x Ksources
S = zeros(Ntrials,Sensors,Tsamples); % Ntrials pairs;
V = zeros(Ntrials,Ksources,Tsamples); % Ntrials Ksources each
Y = zeros(Ntrials,Ksources,Tsamples); % fft(V)
Z = zeros(Ntrials,Sensors,Tsamples); % fft(S)
X = zeros(Ntrials, 2*Tsamples-1);

% Calculations
for j = 1:Ntrials
    for k = 1:Ksources
        V(j,k,:) = sin(2*pi*(k+10*j-10)*t)+0.1*randn(1,Tsamples);
    end
    S(j,:,:) = A*squeeze(V(j,:,:));
    Y(j,:,:) = fft(squeeze(V(j,:,:))')'/length(t);
    Zj = A*squeeze(Y(j,:,:)); % Zj = fft(Sj)
    Z(j,:,:) = Zj; % N x 2 x T
    Zconj = conj(squeeze(Z(j,2,:)));
    X(j,:) = xcorr(squeeze(Z(j,1,:)),Zconj);
end

% Metrics Declarations
C = zeros(Ntrials, 1);
PLV = zeros(Ntrials,1);
PLI = zeros(Ntrials,1);
WPLI = zeros(Ntrials,1);

% Calculations
for j=1:Ntrials
    Cj = mean(X(j,:))/(mean(abs(squeeze(Z(j,1,:)).^2))*...
                       mean(abs(squeeze(Z(j,2,:)).^2)));
    PLVj = mean(exp(1i*angle(X(j,:)))); % usar diff(angle)?
    PLIj = abs(mean(sign(imag(X(j,:)))));
    WPLIj = abs(mean(imag(X(j,:))))/mean(abs(imag(X(j,:)))); 
    
    C(j) = Cj;
    PLV(j) = PLVj;
    PLI(j) = PLIj;
    WPLI(j) = WPLIj;
end




















