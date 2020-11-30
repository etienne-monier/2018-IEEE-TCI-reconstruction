clear all
close all

if ~exist('../3. Final figures/Maps/', 'dir')
   mkdir('../3. Final figures/Maps/')
end
if ~exist('../3. Final figures/Spims/', 'dir')
   mkdir('../3. Final figures/Spims/')
end

%% Load data
addpath('../../algorithms/')
addpath('../../synth_data/')
addpath('../../reproducible_results/')

disp('Loading data ...')

load('synthetic_data.mat','X','Spectra','Elements', 'Maps'),
load('sampling_path')

load('../1. Data2use/Full.mat'), 
load('../1. Data2use/PCA_NLm.mat'), 
load('../1. Data2use/SNN.mat'), 
load('../1. Data2use/SSS.mat'),

%% Structures
FinalMSE = zeros(6,1);

FinalSpectra = cell(6,4);
FianlSAD = zeros(6,4); FianlaSAD = zeros(6,1);

FinalMaps = cell(6,4);
FinalaMapsMSE = zeros(6,1);

FinalSpimDisplay = cell(6,1);

%% Data

[n m l] = size(X);

Data = {X, X2ms, PCA_NLm, XSNN, XSSS};
SampledData = {X, X2ms, X2ms, Y, Y};
ITab = {1:n*m, 1:n*m, 1:n*m, I, I};

%% Normalization

[ Normalisation ] = getNormalisation( PCA_NLm );

%% Init for ground truth

FinalMaps(1,:) = Maps;
FinalSpectra(1,:) = transpose(mat2cell(Spectra',[1 1 1 1],[size(X,3)]));
FinalSpimDisplay{1} = spim2im( X,Normalisation);

%% Processing

disp('Computing data ...')

for n = 1:length(Data)
    
    [M,x,MSE,SAD,aSAD,aMapsMSE,Maps_est] = GetCritera(Data{n},SampledData{n},ITab{n});
    
    FinalMaps(n+1,:) = Maps_est; 
    FinalSpectra(n+1,:) = transpose(mat2cell(M',[1 1 1 1],[size(Data{n},3)])); 
    FinalSAD(n+1,:) = SAD; 
    FinalaSAD(n+1) = aSAD;  
    FinalaMapsMSE(n+1) = aMapsMSE; 
    FinalMSE(n+1) = MSE; 
    FinalSpimDisplay{n+1} = spim2im( Data{n},Normalisation);
end


%% Output figures

disp('Creating figures ...')

% Plot spectra ---
cellYlabel = {'Oracle','Full2ms','PCA_NLm','S2N','3S'};

figure,

% Num of rows and columns
[ fnr,fnc ] = size(FinalSAD);

for ind_x = 2:fnr
    for ind_y = 1:fnc

        subplot(fnr-1,fnc,ind_y + (ind_x-2)*fnc)
        plot([FinalSpectra{ind_x,ind_y}'],'LineWidth',2)
        
        grid on
        
        if (ind_x == 2) title(Elements{ind_y}), end
        if (ind_y == 1) ylabel(cellYlabel{ind_x-1}), end
    
    end
end

% plot maps and save them ---
cellYlabel = {'Actual','Oracle','Full2ms','PCA_NLm','S2N','3S'};

figure,

for ind_x = 1:fnc
    for ind_y = 1:fnr

        subplot(fnc,fnr,ind_y + (ind_x-1)*fnr)
        imshow(uint8(255*FinalMaps{ind_y,ind_x})),
        
        imwrite(uint8(255*FinalMaps{ind_y,ind_x}),strcat('../3. Final figures/Maps/img_',num2str(ind_y),'_',num2str(ind_x),'.png'))
        
        if (ind_y == 1) ylabel(Elements{ind_x},'FontSize',18,'FontWeight','Bold'), end
        if (ind_x == 1) title(cellYlabel{ind_y},'FontSize',18,'FontWeight','Bold'), end
    end
end

% Show spims

cellYlabel = {'Actual','Full','PCA_NLm','SNN','SSS'};
figure,
for ind_y = 1:(fnr-1)
    
    subplot(1,fnr-1,ind_y)
    imshow(FinalSpimDisplay{ind_y+1})
    
    imwrite(FinalSpimDisplay{ind_y+1},strcat('../3. Final figures/Spims/img2_',num2str(ind_y),'.png')) % strcat('../3. Final figures/Spims/',cellYlabel{ind_y},'.png')
    
    title(cellYlabel{ind_y}),
end

%% Latex tabular

disp('Writting latex file ...')

fileID = fopen('../3. Final figures/latexData.tex','w');

CellTitleCell = {'Oracle','Full2ms','PCA_NLm ','S2N	','3S'};

for ind = 1:fnr-1

    fprintf(fileID,'%%\n%%\n');
    
    fprintf(fileID,[CellTitleCell{ind} '\n']);
    fprintf(fileID,'&%.5f\n',FinalMSE(ind+1));
    
    fprintf(fileID,'&%.5f\n',FinalaSAD(ind+1));
    
    fprintf(fileID,'&%.5f\n',FinalaMapsMSE(ind+1));
    fprintf(fileID,'\\\\');
end
fclose(fileID);
















function [ Normalisation ] = getNormalisation( X )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

X = transpose(reshape( X,[size(X,1)*size(X,2) size(X,3)] ));

Bands = [236 346 709];
R = X(Bands(1),:);
G = X(Bands(2),:);
B = X(Bands(3),:);
     

 % third parameter does not exist, so default it to something
 % Normalisation = [ar,br,ag,bg,ab,bb];
 a = @(x) min(x(:));
 b = @(x) max(x(:))-min(x(:));
 Normalisation = [a(R) b(R) a(G) b(G) a(B) b(B)];

end

function [ Im ] = spim2im( X, Normalisation )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Bands = [236 346 709];

NormR = @(x) (x-Normalisation(1))/Normalisation(2);
NormG = @(x) (x-Normalisation(3))/Normalisation(4);
NormB = @(x) (x-Normalisation(5))/Normalisation(6);

Im = cat(3, NormR(X(:,:,Bands(1))) , NormG(X(:,:,Bands(2))) ,NormB(X(:,:,Bands(3))) );
Im(Im>1) = 1;
Im(Im<0) = 0;

end

