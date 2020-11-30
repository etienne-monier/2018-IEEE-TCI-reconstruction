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


%% Create acquired data

disp('Creating Full10ms ...')

[n,m,l] = size(X);

SIGMA = 8.9e-3;

% Noise generation
rng(0);
N = randn(n, m, l);

Y = X + SIGMA*N;

% Sub sampling
Ytemp = transpose(reshape(Y,[n*m l]));
Ytemp(:,~ismember(1:n*m,I)) = 0;
Y  = reshape( Ytemp', [n m l] ); 
clear Ytemp

%% Reconstruction with S2N

disp('Reconstructing image with S2N ...')

rng(1);
InitPoint = randn(l, n * m);

[ XSNN, lambda, mu ] = SNN( Y, I, 'init',InitPoint);

%% Reconstruction with 3S

mu = 1;

disp('Reconstructing image with 3S ...')

XSSS = SSS( Y, I, mu, 'init', InitPoint ,'verbose','n' );

%% Saving data

disp('Saving data ...')

if ~exist('../1. Data2use/', 'dir')
   mkdir('../1. Data2use/')
end

save '../1. Data2use/SSS.mat' XSSS
save '../1. Data2use/SNN.mat' XSNN lambda mu