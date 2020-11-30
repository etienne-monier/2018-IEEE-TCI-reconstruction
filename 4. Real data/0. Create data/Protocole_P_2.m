clear all
clc

addpath('../../algorithms/')
addpath('../../real_data')

%% Load data

disp('Loading data ...')

load('Full10ms.mat')

[n,m,l] = size(Full10ms);

rng(1);
InitPoint = randn(l, n * m);

load('path.mat')
load('seed.mat')

%% Setup

disp('Creating Partial ...')

PIX_RATIO = 0.2;

ns = floor(PIX_RATIO*n*m);
I = sub2ind([n,m], yScan(1:ns), xScan(1:ns));
A = zeros([l n*m]); A(:,I) = 1; A = reshape(A',[n m l]);
Y = Full10ms.*A;

%% Reconstruct with 3S

disp('Reconstructing image with 3S ...')

lambda = 0.001;

[ XSSS ] = SSS( Y, I, lambda, 'init',InitPoint,'verbose','n');

%% Reconstruct with S2N

[ XSNN, lambda, mu ] = SNN( Y, I, 'init',InitPoint);

%% Save results

disp('Saving ...')

if ~exist('../1. Data2use/', 'dir')
   mkdir('../1. Data2use/')
end

save '../1. Data2use/XSNN.mat' XSNN lambda mu
save '../1. Data2use/XSSS.mat' XSSS
save '../1. Data2use/partial_path.mat' I
