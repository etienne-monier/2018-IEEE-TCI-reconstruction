clear all
clc

addpath('../algorithms/')
addpath('../synth_data/')
addpath('../reproducible_results/')

%% Load data

load('synthetic_data.mat')
load('sampling_path')

[n,m,l] = size(X);

%% Reproducible results

% Noise matrix
rng(0);
Noise = randn(n,m,l);

% Init point
rng(1);
InitPoint = randn(l, n * m);

%% Parameters

SNR = 25.1;
PIX_RATIO = 0.2;

%% Observation matrix creation

% Noise standard deviation
Pow = mean(X(:).^2);
SIGMA = sqrt(Pow*10^(-SNR/10));

% Observation matrix
Y = X + SIGMA*Noise;

%% Sub sampling

% Sampling path
% I comes from sampling_path.mat

Ytemp = transpose(reshape(Y,[n*m l])); Ytemp(:,~ismember(1:n*m,I)) = 0;
Y  = reshape( Ytemp',[n m l] ); 
clear Ytemp

%% Getting the parameter matrix

MuTab = 10.^[-5:1:5];

%% Reconstruction

ind = 1;
for mu = MuTab

    fprintf('mu = %.5f (it %d/%d)\n',mu,ind,length(MuTab));
    
    Xhat = SSS( Y, I, mu, 'init', InitPoint ,'verbose','n' );
    
    NMSE(ind) = norm(Xhat(:) - X(:), 'fro')^2 / (norm(X(:), 'fro')^2);
    
    ind = ind + 1;
end

%% Save results

save 'parameter_search_SSS.mat' MuTab NMSE;

%% Display results

figure,
loglog(MuTab, NMSE)

xlabel('\mu'), ylabel('NMSE')
title('3S reconstruction performances w.r.t. \mu')
