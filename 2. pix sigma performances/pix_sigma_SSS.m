clear all
clc

addpath('../algorithms/')
addpath('../synth_data/')

%% Data

load('synthetic_data.mat', 'X')

[n,m,l] = size(X);
Pow = mean(X(:).^2);


%% Reproducible results

% Full sampling path
rng(1);
I0 = randperm(n * m);

% Init point
rng(1);
InitPoint = randn(l, n * m);


%% Parameter tabs

Nrel = 10;                  % # realizations
SnrTab = 19:2:29;           % Snr tab
PixTab = [0.2 0.3 0.4];     % pixel ratio tab

itTab = 1:Nrel;             % iteration tab

%% Process

% Data stuctures
parameters = combvec(itTab,SnrTab,PixTab);
NMSETemp = zeros(1,size(parameters,2));
MuTemp = zeros(1,size(parameters,2));

% Create pool
p = gcp(); % If no pool, do not create new one.
disp(['Pool of size ' num2str(p.NumWorkers) ' used'])

parfor ind = 1:size(parameters,2)
    
    % Print the status and save it to a file
    fprintf('ind = %d\n',ind)
    fileID = fopen('status.txt','a');
    fprintf(fileID,'ind = %d\n',ind);
    fclose(fileID);
    
    % Parameters to be used
    it = parameters(1,ind);
    SNR = parameters(2,ind);
    pix_r = parameters(3,ind);
    
    % Data init
    Noise = randn(n,m,l);
    SIGMA = sqrt(Pow*10^(-SNR/10));
    Y = X + SIGMA*Noise;
    I = sort(I0(1:floor(pix_r*n*m)));
    Ytemp = transpose(reshape(Y,[n*m l])); Ytemp(:,~ismember(1:n*m,I)) = 0;
    Y = reshape( Ytemp',[n m l] );
    
    % Search of the optimum parameter for 3S
    MuTab = 10.^[-1:0.1:0];

    MSE = zeros(1,length(MuTab));
    ind0 = 1;
    for mu = MuTab
        Xhat = SSS( Y, I, mu, 'init', InitPoint ,'verbose','n' );
        MSE(ind0) = norm(Xhat(:) - X(:), 'fro')^2 / (norm(X(:),'fro')^2);
        ind0 = ind0 + 1;
    end

    % Final optimal parameter and NMSE
    MuTemp(ind) = MuTab(find(MSE == min(MSE)));
    NMSETemp(ind) = min(MSE);
    
end

% Reshape output
NMSE = permute(reshape(NMSETemp,[Nrel,length(SnrTab),length(PixTab)]), [3 2 1]);
Mu = permute(reshape(MuTemp,[Nrel,length(SnrTab),length(PixTab)]) , [3 2 1]);

%% Save results

if ~exist('output', 'dir')
   mkdir('output')
end

save 'output/pix_sigma_SSS_results.mat' SnrTab PixTab NMSE Mu

%% Display

Err = std(NMSE, 0, 3);

Color = {'blue','red','green'};

figure,
hold on
for ind = 1:length(PixTab)
    e = errorbar(SnrTab, mean(NMSE(ind,:,:),3) ,Err(ind,:));
    e.Marker = '*';
    e.MarkerSize = 10;
    e.Color = Color{ind};
    e.CapSize = 15;
    grid on
    
end

xlabel('SNR (dB)'), ylabel('NMSE'),
legend('r = 0.2', 'r = 0.3', 'r = 0.4'),
title('Performance of 3S in term of NMSE w.r.t. the pixel ratio r and noise level')
