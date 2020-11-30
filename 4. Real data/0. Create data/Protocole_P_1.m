clear all
clc

addpath('../../algorithms/')
addpath('../../real_data')

%% Load data

disp('Loading data ...')

load('Full2ms.mat')
[n,m,l] = size(Full2ms);

%% Setup

disp('Denoising ...')
[ Xout ] = PCAdenoise( Full2ms, 20 );

[ SIGMA ] = madNoiseEst( transpose(reshape(Xout,[n*m l])) );

[PCA_NLm] = hyper_denoising_nlmeans(Xout,[],[1],5*SIGMA,SIGMA);

disp('Saving data ...')

if ~exist('../1. Data2use/', 'dir')
   mkdir('../1. Data2use/')
end

save '../1. Data2use/PCA_NLm.mat' PCA_NLm
