%% QuintoMes - Granger Causality (prediction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

% setup
clear, clc, clf
% load('.\Data.mat\SpkBuz.mat')
% load('.\Data.mat\thetaBuz.mat')
load('.\Data.mat\LFPBuz.mat')
LFP = LFP';

% estimate autocovariance sequence from time series data
[A2,Sig,E2]= tsdata_to_var([LFP],15); %q=15

% Returns autocovariance sequence G defined as
%%%%%%%      L_k = cov(X_t, X_(t-k))      %%%%%%%
% for a VAR model with coefficients A and
% (positive-definite) residual covariance matrix
% SIG, by "reverse-solving" the Yule-Walker
% equations
%%%%% L_k = sum{P}_{l=1}(A_l*L_(k-l) + del_k0*SIG %%%%%
% (where  = SIG). The algorithm solves the
% associated 1-lag problem - a discrete-time
% Lyapunov equation - and then calculates higher
% lags recursively [1].    
[G,info] = var_to_autocov(A2,Sig);

% Displays information (errors, warnings,
% diagnostics, etc.) reported by the routine
% var_to_autocov; see that routine for details  
acerr = var_acinfo(info, 0) % abort_on_error = false

% Returns (column) vector freqs of fres+1 equally
% spaced frequencies on [0,fs/2], where fs is a
% sample rate (so fs/2 is the Nyqvist frequency).
% If a sample rate is not supplied (default),
% frequencies are returned over the normalised
% range [0 pi].     
freqs = sfreqs(1000,srate); %freq_res=1000

% Returns the matrix f of pairwise-conditional
% frequency-domain (spectral) MVGCs 
%%%%% f_ij(lam) = f_(X_j-->X_i|X_[ij])(lam) %%%%%
% (where [ij] denotes omission of the ij-th
% variables) between all pairs of variables
% represented in G, for a stationary VAR process
% with autocovariance sequence G. The first index
% i of f is the target (causee) variable, the
% second j the source (causal) variable and the
% third indexes the frequency      
[F_spect,fres] = autocov_to_spwcgc(G,1000,[]); %freq_res=1000

ts = squeeze(F_spect);

plot(1:length(g),g)


%% Real Data
% Parameters

clear, clf, clc
ntrials   = 10;     % number of trials
nobs      = 1000;   % number of observations per trial

regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = 'F';    % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDRD'; % multiple hypothesis test correction (see routine 'significance')

fs        = 1250;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

seed      = 0;      % random seed (0 for unseeded)

%% Data
load('EE.210.mat')
pfc = data(11:12,12500:10:75000);
ca = data(71:72,12500:10:75000);
% plot(pfc'), hold on
% plot(ca'), hold off

% fazer zscore e downsampling 
X = [pfc; ca];
% X = data(1:96,12500:75000);

%% Model order estimation (<mvgc_schema.html#3 |A2|>)

% Calculate information criteria up to specified maximum model order.

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria.

figure(1); clf;
plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
title('Model order estimation');

AT = tsdata_to_var(X,15,regmode)
amo = size(AT,3); % actual model order

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);
fprintf('actual model order     = %d\n',amo);

% Select model order.

if     strcmpi(morder,'actual')
    morder = amo;
    fprintf('\nusing actual model order = %d\n',morder);
elseif strcmpi(morder,'AIC')
    morder = moAIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
elseif strcmpi(morder,'BIC')
    morder = moBIC;
    fprintf('\nusing BIC best model order = %d\n',morder);
else
    fprintf('\nusing specified model order = %d\n',morder);
end

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)

% Estimate VAR model of selected order from data.

ptic('\n*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(X,morder,regmode);
assert(~isbad(A),'VAR estimation failed - bailing out');
ptoc;

% Report information on the estimated VAR, and check for errors.
%
% _IMPORTANT:_ We check the VAR model for stability and symmetric
% positive-definite residuals covariance matrix. _THIS CHECK SHOULD ALWAYS BE
% PERFORMED!_ - subsequent routines may fail if there are errors here. If there
% are problems with the data (e.g. non-stationarity, colinearity, etc.) there's
% also a good chance they'll show up at this point - and the diagnostics may
% supply useful information as to what went wrong.

info = var_info(A,SIG);
assert(~info.error,'VAR error(s) found - bailing out');

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities from VAR model parameters
% by state-space method [4]. The VAR model is transformed into an equivalent state-
% space model for computation. Also return p-values for specified test (F-test or
% likelihood-ratio test; this is optional - if p-values are not required, then it
% is not necessary to supply time series |X|, regression mode |regmode|, or test
% specification |tstat|).

ptic('*** var_to_pwcgc... ');
[F,pval] = var_to_pwcgc(A,SIG,X,regmode,tstat);
ptoc;

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed - bailing out');

% Significance-test p-values, correcting for multiple hypotheses.

sig = significance(pval,alpha,mhtc);

% Plot time-domain causal graph, p-values and significance.

figure(2); clf;
sgtitlex('Pairwise-conditional Granger causality - time domain');
subplot(1,3,1);
plot_pw(F);
title('Pairwise-conditional GC');
subplot(1,3,2);
plot_pw(pval);
title(['p-values (' tstat '-test)']);
subplot(1,3,3);
plot_pw(sig);
title(['Significant at \alpha = ' num2str(alpha)]);

%% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)

% If not specified, we set the frequency resolution to something sensible. Warn if
% resolution is very large, as this may lead to excessively long computation times,
% and/or out-of-memory issues.

if isempty(fres)
    fres = 2^nextpow2(info.acdec); % based on autocorrelation decay; alternatively, you could try fres = 2^nextpow2(nobs);
	fprintf('\nfrequency resolution auto-calculated as %d (increments ~ %.2gHz)\n',fres,fs/2/fres);
end
if fres > 20000 % adjust to taste
	fprintf(2,'\nWARNING: large frequency resolution = %d - may cause computation time/memory usage problems\nAre you sure you wish to continue [y/n]? ',fres);
	istr = input(' ','s'); if isempty(istr) || ~strcmpi(istr,'y'); fprintf(2,'Aborting...\n'); return; end
end

% Calculate spectral pairwise-conditional causalities at given frequency
% resolution by state-space method.

ptic('\n*** var_to_spwcgc... ');
f = var_to_spwcgc(A,SIG,fres);
assert(~isbad(f,false),'spectral GC calculation failed - bailing out');
ptoc;

% Plot spectral causal graph.

figure(3); clf;
sgtitlex('Pairwise-conditional Granger causality - frequency domain');
plot_spw(f,fs);

%% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)

% Check that spectral causalities average (integrate) to time-domain
% causalities, as they should according to theory.

fprintf('\nfrequency-domain GC integration check... ');
Fint = smvgc_to_mvgc(f); % integrate spectral MVGCs
amax = maxabs(F+Fint)/2;
if amax < 1e-5; amax = 1; end % in case all GCs very small
mre = maxabs(F-Fint)/amax;
if mre < 1e-5
    fprintf('OK (maximum relative error ~ %.0e)\n',mre);
else
    fprintf(2,'WARNING: high maximum relative error ~ %.0e\n',mre);
end



















