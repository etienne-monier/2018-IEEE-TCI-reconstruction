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

Ytemp = transpose(reshape(Y,[n*m l])); Ytemp(:,~ismember(1:n*m,I)) = 0;
Y  = reshape( Ytemp',[n m l] ); 
clear Ytemp

%% Getting the parameter matrix

% Parameter tab
LambdaTab = 10.^[-4:0.25:2];
MuTab = 10.^[-5:0.25:3];

% Parameter Tab
parameters = combvec(LambdaTab, MuTab);

Nparam = length(parameters);
Nlambda = length(LambdaTab);
Nmu = length(MuTab);

%% Reconstruction

for ind = 1:Nparam

    fprintf('lambda = %.5f and mu = %.5f (ind %d/%d)\n', parameters(1, ind), parameters(2, ind), ind, Nparam);

    % Performing reconstruction
    Xhat = SNN(Y, I, 'init', InitPoint , 'verbose', 'y', 'parameters', parameters(:, ind));
    
    % Storing performance.    
    NMSE(ind) = norm(Xhat(:)-X(:), 'fro')^2 / (norm(X(:), 'fro')^2);

end

NMSE = reshape(NMSE, [Nlambda Nmu]);

%% Save results

save 'parameter_search_SNN.mat' LambdaTab MuTab NMSE;


%% Display results

% Getting colormap
T = real2rgb(log(NMSE), 'parula');

% Displaying figure
figure,
surf(LambdaTab, MuTab, log10(NMSE'), T),

xlabel('lambda'),ylabel('mu'),zlabel('MSE'),
set(gca, 'XScale', 'log', 'YScale', 'log'),
view([0 0 90]),

c = colorbar;
c.Label.String = 'log(MSE)';


% hold on,

% h = 0;
% caxis([min(log10(temp(:))) max(log10(temp(:)))]);

% % First plot : Lambda alone
% h1 = plot3([EmpiricChoice.LambdaAlone EmpiricChoice.LambdaAlone],[MuTab(1) MuTab(end)],[h h]);
% h1.LineWidth = 6;
% h1.LineStyle = '--';
% h1.Color = 'red';


% % Second plot : Mu alone
% h2 = plot3([LambdaTab(1) LambdaTab(end)],[EmpiricChoice.MuAlone EmpiricChoice.MuAlone],[h h]);
% h2.LineWidth = 6;
% h2.LineStyle = ':';
% h2.Color = [204, 51, 255]/255;

% % Third plot : both Lambda and Mu
% h3 = scatter3([EmpiricChoice.FinalLambda],[EmpiricChoice.FinalMu],[h],150,'filled');
% h3.MarkerEdgeColor = 'black';
% h3.MarkerFaceColor = 'white';

% % Minimum on grid
% [m, indices] = minn(temp);
% h4 = scatter3([LambdaTab(indices(2))],[MuTab(indices(1))],[h],150,'filled');
% h4.MarkerEdgeColor = 'white';
% h4.MarkerFaceColor = 'black';

% legend([h1 h2 h3 h4],{'\lambda_*','\mu_*','Chosen parameters','Optimal parameters on grid'},'Location','northwest','FontSize',18)


