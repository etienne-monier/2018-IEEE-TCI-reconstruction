clear all
clc

addpath('../../real_data')

%% Load data

disp('Loading data ...')

load('Full10ms.mat')
[n, m, l] = size(Full10ms);

load('path.mat')

PIX_RATIO = 0.2;

%% Setup

disp('Creating Partial ...')

ns = floor(PIX_RATIO*n*m);
I = sub2ind([n,m], yScan(1:ns), xScan(1:ns));
A = zeros([l n*m]); A(:,I) = 1; A = reshape(A',[n m l]);

Partial = Full10ms.*A;

%% Save

disp('Saving ...')

if ~exist('../1. Data2use/', 'dir')
   mkdir('../1. Data2use/')
end

save '../1. Data2use/Partial.mat' Partial
