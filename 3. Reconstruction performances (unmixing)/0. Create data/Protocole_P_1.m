clear all
clc

addpath('../../algorithms/')
addpath('../../synth_data/')
addpath('../../reproducible_results/')

%% Load data
disp('Loading data ...')

load('synthetic_data.mat', 'X')

% Reproducible results
load('sampling_path')

% Synthetic data
[n,m,l] = size(X);

%% Full 2ms

disp('Creating Full2ms ...')

% Add noise (2ms)
SNR = 19.22;
Pow = mean(X(:).^2);
SIGMA = sqrt(Pow*10^(-SNR/10));

% Noise generation
rng(0);
N = randn(n, m, l);

X2ms = X + SIGMA*N;

%% PCA + NLm

disp('Creating PCA_NLm ...')
Xout = PCAdenoise( X2ms, 20 );

SIGMA = madNoiseEst( transpose(reshape(Xout,[n*m l])) );

PCA_NLm = hyper_denoising_nlmeans(Xout,[],[1],0.2*SIGMA,SIGMA);

%% 3S

lambda = 0.001;
disp('3S denoising ...')

rng(1);
InitPoint = randn(l, n * m);

[ XSSS ] = SSS( X2ms, 1:m*n, lambda, 'init',InitPoint,'verbose','n');

%% S2N

disp('S2N denoising ...')
[ XSNN_den, lambda, mu ] = SNN( X2ms, 1:m*n, 'init',InitPoint);


%% Save data

disp('Saving data ...')

if ~exist('../1. Data2use/', 'dir')
   mkdir('../1. Data2use/')
end

save '../1. Data2use/Full2ms.mat' X2ms
save '../1. Data2use/PCA_NLm.mat' PCA_NLm
save '../1. Data2use/3S_denoising.mat' XSSS_den
save '../1. Data2use/S2N_denoising.mat' XSNN_den lambda mu
