%% Primeiro estudo de coerencia
% 1) Efeito da qtd de win na coerencia com muitas 
% wins e com poucas wins
% 2) Efeito do tamanho da winindows
% 3) Efeito da sobreposição de windows ****
% 4) verificar a significancia estatistica de cada 
% um desses efeitos.
% 5) existe um trade off entre o tamanho e a quantidade
% de janelas janelas maiores tornan o transformada de 
% fourier mais precisa (dado estacionario), mas 
% diminuindo a quantidade de janelas o vies da
% quantidade de janelas aumenta.


% em um experimento o animal performa durante 
% um determinado tempo. Pode ser uma tarefa 
% longa realizada poucas vezes ou pode ser uma
% tarefa curta realizada varias vezes. Posso 
% utilizar varios animais ou poucos.

% uma terefa longa dura 25min
% uma tarefa curta dura 30s
% muitos animais são 10 e poucos são 3

%% CELL(0) Testando 2 e 3
% % usando a funcao built-in do matlab **CODIGO DA AULA 6**
% mscohere % ms = magnitude squared

%%% rodar muitas vezes sem limpar a figura e alterando o 
%%% Tmax para estudar o efeito do tamanho da amostra.


%setup
clear, clc, clf
srate = 1000;
dt = 1/srate;
% alterar o tamanho da amostra muda a coerencia
% pois altarea a quantidade de vetores sobre os 
% quais eu estou computando a media

figure(4)
subplot(211)
Tmax = 2;
cont = 0;
while (Tmax <= 32)
Tmax = Tmax*2
t = dt:dt:Tmax;
cont = cont+1;

% o valor randn add não faz diferença pois ele não 
% esta mudando over time
phi = -deg2rad(90)+0*randn;

% se retirar o randn de X e Y teremos coerencia 1 em 
% todas as frequencias isso é o que está acontecendo no programa
X = 3*sin(2*pi*10*t)+0.3*randn(size(t));
Y = 0.5*sin(2*pi*10*t+phi)+0.3*randn(size(t));

% aumentar o tamanho da janela gera mais precisão na freq
% mas gera mais ruido nas outras.
windowlength = 2*srate;
overlap = 0;

nfft = 2^16;
[Cxy F] = mscohere(X,Y,windowlength,overlap,nfft,srate);
hold on
plot(F,Cxy, 'linew', cont/2)
xlabel('Frequency (Hz)')
ylabel('Coherence')
xlim([0 20])
ylim([0 1])
end

legend('4', '8', '16', '32', '64')
title('Effect of window quantity')
hold off

% %

%setup
clear
srate = 1000;
dt = 1/srate;

subplot(212)
Tmax = 64;
cont = 32;
contt = 0;
while (Tmax/cont < 64)
t = dt:dt:Tmax;
cont = cont/2
contt = contt + 1;

phi = -deg2rad(90)+0*randn;
X = 3*sin(2*pi*10*t)+0.3*randn(size(t));
Y = 0.5*sin(2*pi*10*t+phi)+0.3*randn(size(t));

windowlength = cont*srate;
overlap = 0;

nfft = 2^16;
[Cxy F] = mscohere(X,Y,windowlength,overlap,nfft,srate);
hold on
plot(F,Cxy, 'linew', contt/2)
xlabel('Frequency (Hz)')
ylabel('Coherence')
xlim([0 20])
ylim([0 1])
end

legend('16', '8', '4', '2', '1')
title('Effect of window size')
hold off

%% Cell(1) Média vet de vet unitarios como previwe da coh
% cria N vetores unitarios e computa a coherencia
% se eu pegar cada vez mais janelas eu não estou 
% vendo o trad off das janelas serem menores
% assumindo que o tamanho do sinal é fixo.
% tudo que estou vendo é que mais win dá mais acurácia

clc % limpa comm win
clf % limpa fig
clear % limpa mem

Nvetor = 100; % numero de vetores
Ni = 10000; % numero de simulações
NAlliAll = zeros(Ni, Nvetor);
% usar uma seed para os dados randn

tic
for nvetor = 1:Nvetor
    nvetor;
    for i = 1:Ni
        clear theta Vall 
        theta = pi/2 + 0.3*randn*(1:nvetor); 
        Vall = exp(1i*theta);
        
        Cxy = abs(mean(Vall));
        NAlliAll(i,nvetor) = Cxy;
    end
end
toc % Elapsed time is 36.704769 seconds.

V = var(NAlliAll);
figure(1)
plot(1:Nvetor,V)
xlabel('# windows')
ylabel('Varience')

NAlliAllOfDirectVector03randAmp1000simulacao75vetores = NAlliAll;

%% Cell(2) As Cell(1) but with Imag(coh)

clc, clf, clear

Nvetor = 75; % numero de vetores
Ni = 1000; % numero de simulações
NAlliAll = zeros(Ni, Nvetor);
% usar uma seed para os dados randn

tic
for nvetor = 1:Nvetor
    nvetor;
    for i = 1:Ni
        clear theta Vall 
        theta = pi/2 + 1*randn*(1:nvetor); 
        Vall = exp(1i*theta);
        
        Cxy = abs(imag(mean(Vall)));
        NAlliAll(i,nvetor) = Cxy;
    end
end
toc % Elapsed time is 3.135485 seconds.

V = var(NAlliAll);
figure(1)
plot(1:Nvetor,V)
xlabel('# windows')
ylabel('Varience')

NAlliAllOfDirectVectImag1randAmp1000simulacao75vetores = NAlliAll;

%% Cell(3) As Cell(1) but changing noise amp
% ver se existe uma constante / proporção clara 
% entre o ruido no angulo e a curva de decaimento

clc % limpa comm win
clf % limpa fig
clear % limpa mem

Nvetor = 100; % numero de vetores
Ni = 1000; % numero de simulações
Amp = 0.01:0.02:2; % 100 amplitudes de ruido
NAlliAll = zeros(Ni, Nvetor, length(Amp));
% usar uma seed para os dados randn

tic
cont = 0;
for amp = Amp
    cont = cont+1
    for nvetor = 1:Nvetor
        for i = 1:Ni
            clear theta Vall 
            theta = pi/2 + amp*randn*(1:nvetor); 
            Vall = exp(1i*theta);

            Cxy = abs(mean(Vall));
            NAlliAll(i,nvetor,cont) = Cxy;
        end
    end
end
toc % Elapsed time is 732.839375 seconds.
V = squeeze(var(NAlliAll));

figure(1)
for cont = 2:10:100
plot(1:100,V(cont,:), 'linew', 2)
hold on
end
legend('.01', '0.2', '0.4', '0.6', '0.8', ...
    '1.0', '1.2', '1.4', '1.6', '1.8')
hold off
xlabel('# windows')
ylabel('Varience')
title('Noise Influence in Coherence Varience per Number os Windows in Data')

%% Cell(3) Smothing: moving average by convolution

clear, clc, clf
%loading V from previews lest code block
load('PrimeiroMesCell3.mat')
% ordem = #points no kernel
ordem = 5;
% kernel for moving average
kernel = ones(1,ordem)/ordem;
Convol = zeros(100, 100);
slope = zeros(1,100);
lo=1;
hi=4;
figure(1)
subplot(211)
    for i = 1:100
        Convol(i,:) = conv(V(i,:),kernel,'same');
        slope(i) = (Convol(i,lo)-Convol(i,hi))/(hi-lo+1);
        plot(1:100,Convol(i,:)), hold on
        plot(lo:hi,Convol(i,lo:hi),'yo')
    end
    xlim([lo-5 hi+5])
    hold off
subplot(212)
    plot(Amp,slope)
    n = 3; % grau do polinomio
    [p, S] = polyfit(Amp,slope,n);
    Y1 = polyval(p,Amp);
    hold on
    plot(Amp,Y1,'g-','linew',3)
    hold off

figure(2)
for cont = 2:10:100
plot(1:100,Convol(cont,:), 'linew', 2)
hold on
end
legend('.01', '0.2', '0.4', '0.6', '0.8', ...
    '1.0', '1.2', '1.4', '1.6', '1.8')
hold off
xlim([0 95])
xlabel('# windows')
ylabel('Varience')
title('Noise Influence in Coherence Varience per Number os Windows in Data')


%% Cell(4) As Cell(1) but changing #Trials
% Justification of the #Trials by varience stability

clc % limpa comm win
clf % limpa fig
clear % limpa mem

Nvetor = 100; % #vetores
Ni = 10000; % #simulações
V = zeros(4,Nvetor);
Ni = Ni*1;

figure(1) % < 1min
for diff = 1:3 % not reliable with 4
    newNi = Ni/(10^diff); % redução no #simulações
    NAlliAll = zeros(newNi, Nvetor);
    for nvetor = 1:Nvetor
        for i = 1:newNi
            clear theta Vall
            theta = pi/2 + 0.3*randn*(1:nvetor);
            Vall = exp(1i*theta);

            Cxy = abs(mean(Vall));
            NAlliAll(i,nvetor) = Cxy;
        end
    end
    V(diff,:) = var(NAlliAll);
end
plot(1:Nvetor,V(:,:))
xlabel('# windows')
ylabel('Varience')
legend('10000','1000','100')

% X = 1:Nvetor;
% Y = V(1,:);
% p = polyfit(X,Y,8);
% Y = polyval(p,X);
% % 
% hold on
% plot(X,Y,'g-','linew',2)
% hold off

%% Cell(5) Using senosoidal signal
% computa um ponto coh over (fake) windows na freq ff (cell 1) 
% para todas as freq em freqvector.

% como não obtive diferenças na divergencia da coh
% para o tamanho da win em cada teste posso dizer que ...?
% para responder preciso avaliar não apenas a divergencia 
% entre os valores de coh obtidos, mas sim seus valores?
% mas se eles não divergem então os valores entre coh 
% de um mesmo teste não importa. preciso comparar o valor
% da coerencia de um teste com mais janelas com outro com menos


clear, clf, clc
srate = 1000; f = 15; ff = 25; %[15 25]
dt = 1/srate; Tmax = 4; t = dt:dt:Tmax;

Nvetor = 75; % #vetores
Ni = 100; % #simulações
NAlliAll = zeros(Ni, Nvetor);
W = hamming(length(t))'; %vai pra fora
K = exp(-1i*2*pi*ff*t);
tic
for nvetor = 1:Nvetor
    nvetor;
    for i = 1:Ni
        clear nFXYAll
        for nwindow = 1:nvetor %fake signal de N*1seg.
            phi = -deg2rad(90)+0.3*randn;
            X = 3*sin(2*pi*f*t)+0.3*randn(size(t));
            Y = 0.5*sin(2*pi*f*t+phi)+0.3*randn(size(t));
            FX = mean((W.*X).*K); 
            FY = mean((W.*Y).*K);
            nFX = FX/abs(FX);
            nFY = FY/abs(FY);
            nFXY = nFX*conj(nFY);
            Vall(nwindow) = nFXY;
        end        
        Cxy = abs(mean(Vall));
        NAlliAll(i,nvetor) = Cxy;
    end
end
toc % Elapsed time is 308.218937 seconds. Tmax = 1
    % Elapsed time is 539.868770 seconds. Tmax = 2
    % Elapsed time is 954.092664 seconds. Tmax = 4
%%%%%%%%%%%%%%%%%%%%%% aumentar para 100s %%%%%%%%%%%%%%%%%%%%%%%
V = var(NAlliAll);
plot(1:Nvetor,V)
xlabel('# windows')
ylabel('Varience')

NAlliAll15F25FF2segwin1000simulacao75vetores = NAlliAll;


%% Cell(6) As Cell(5) but with fft for speed

% setup
clear, clf, clc
srate = 1000; f = 15; nyquist = 1/srate;
dt = 1/srate; Tmax = 2; t = dt:dt:Tmax;
N = length(t);
F = linspace(0,nyquist,floor(N/2)+1);

Nvetor = 75; % #vetores
Ni = 100; % #simulações

%%% Nolte com many trials %%%
tic
Fx = zeros(Nvetor,N);
Fy = zeros(Nvetor,N);
NAlliAll = zeros(Ni,Nvetor,N);
for nvetor = 1:Nvetor
    nvetor;
    for i = 1:Ni
        for nwin = 1:nvetor
            phaseLag = deg2rad(pi/2)+0.3*randn;
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
        % preencho coluna a coluna colocando um vetor como elemento
        NAlliAll(i,nvetor,:) = abs(C);% diff de (5) aqui guardo um vet
    end
end
toc % Elapsed time is 36.717146 seconds. Tmax = 1s
    % Elapsed time is 85.532558 seconds. Tmax = 2s 
    % Elapsed time is 166.517847 seconds. Tmax = 4s 
V = zeros(Nvetor,length(F));
for ff = 1:length(F)
    V(:,ff) = var(NAlliAll(:,:,ff));
end
for ff = 2:length(F)-1 % a ff 1 e 2001 são diferentes das outras
    plot(1:Nvetor,V(:,ff))
    hold on
end
hold off
xlabel('# windows')
ylabel('Varience')

%% Cell(7) As Cell(6) but with Real Data
% computa um ponto coh over windows na freq ff 
% para todas as freq em freqvector.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% NÃO TERMINADO %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1,2,5,10,25,50,75 nwi

% setup
clear, clf, clc, load('EE.210.mat')
srate = sampleRate; nyquist = 1/srate; clear sampleRate
dt = 1/srate; Tmax = length(data)/srate; t = dt:dt:Tmax;

N = length(t);
F = linspace(0,nyquist,floor(N/2)+1);

Nvetor = 75; % #vetores

%%% Nolte com many trials %%%
tic
Fx = zeros(Nvetor,N);
Fy = zeros(Nvetor,N);
NAlliAll = zeros(Nvetor,N);
for nvetor = 1:Nvetor
    nvetor;
    for nwin = 1:nvetor
        phaseLag = deg2rad(pi/2)+0.3*randn;
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
    % preencho coluna a coluna colocando um vetor como elemento
    NAlliAll(i,nvetor,:) = abs(C);% diff de (5) aqui guardo um vet
end
toc 
V = zeros(Nvetor,length(F));
for ff = 1:length(F)
    V(:,ff) = var(NAlliAll(:,:,ff));
end
for ff = 2:length(F)-1 % a ff 1 e 2001 são diferentes das outras
    plot(1:Nvetor,V(:,ff))
    hold on
end
hold off
xlabel('# windows')
ylabel('Varience')

%% Cell(8) Variancia da coh do sinal todo considerando diferentes qtd de win
% cada ponto da matriz NAlliAll possui o valor da coh do sinal
% O tipo de analise esta dividido nas linhas da matriz
% cada linha possui um numero diferente de janelas analisadas
% foi repetida essa analise uma vez para cada coluna

clc % limpa comm win
clf % limpa fig
clear % limpa mem

load('NAlliAll15F15FF1segwin1000simulacao75vetores.mat')
load('NAlliAll15F25FF1segwin1000simulacao75vetores.mat')
load('NAlliAll15F25FF2segwin1000simulacao75vetores.mat')
load('NAlliAll15F25FF4segwin1000simulacao75vetores.mat')
load('NAlliAllOfDirectVector1randAmp1000simulacao75vetores.mat')
load('NAlliAllOfDirectVector03randAmp1000simulacao75vetores.mat')

% Variancia das analises
V = zeros(6,75);
V(1, :) = var(NAlliAll15F15FF1segwin1000simulacao75vetores);
V(2, :) = var(NAlliAll15F25FF1segwin1000simulacao75vetores);
V(3, :) = var(NAlliAll15F25FF2segwin1000simulacao75vetores);
V(4, :) = var(NAlliAll15F25FF4segwin1000simulacao75vetores);
V(5, :) = var(NAlliAllOfDirectVector1randAmp1000simulacao75vetores);
V(6, :) = var(NAlliAllOfDirectVector03randAmp1000simulacao75vetores);

figure(1)
for cont = 1:6
plot(1:75,V(cont,:))
hold on
end
xlabel('# windows')
ylabel('Varience')
legend('15F15FF1segwin', '15F25FF1segwin', '15F25FF2segwin', ...
    '15F25FF4segwin', 'DirectVector1rand', 'DirectVector03rand')
hold off

%% Cell(9) visualização dos valores da coh entre diff win length

M15F25FF1seg = NAlliAll15F25FF1segwin1000simulacao75vetores;
M15F25FF2seg = NAlliAll15F25FF2segwin1000simulacao75vetores; 
M15F25FF4seg = NAlliAll15F25FF4segwin1000simulacao75vetores;

figure(1)
subplot(2,2,1)
    plot(M15F25FF1seg(20,10:75),M15F25FF2seg(20,10:75),'ko')
    x = M15F25FF1seg(20,10:75);
    y = M15F25FF2seg(20,10:75);
    n = 1; % grau do polinomio
    p = polyfit(x,y,n);
    X = 0:0.01:0.6;
    Y1 = polyval(p,X);
    hold on
    plot(X,Y1,'g-','linew',3)
    hold off
    xlim([0 .6])
    ylim(xlim)
    axis square
    
subplot(2,2,2)
plot(M15F25FF1seg(20,10:75),M15F25FF4seg(20,10:75),'ko')
    x = M15F25FF1seg(20,10:75);
    y = M15F25FF4seg(20,10:75);
    p = polyfit(x,y,n);
    X = 0:0.01:0.6;
    Y2 = polyval(p,X);
    hold on
    plot(X,Y2,'g-','linew',3)
    hold off
    xlim([0 .6])
    ylim(xlim)
    axis square
    
subplot(2,2,3)
    plot(M15F25FF2seg(20,10:75),M15F25FF4seg(20,10:75),'ko')
    x = M15F25FF2seg(20,10:75);
    y = M15F25FF4seg(20,10:75);
    p = polyfit(x,y,n);
    X = 0:0.01:0.6;
    Y = polyval(p,X);
    hold on
    plot(X,Y,'g-','linew',3)
    hold off
    xlim([0 .6])
    ylim(xlim)
    axis square

    
    %verificar melhor se ha diferença na coh
    % de sinais de janelas de tamanhos diferentes
subplot(2,2,4)
diff = Y1-Y2;
plot(1:61, diff)






