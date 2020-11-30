clear all
close all

if ~exist('../3. Final figures/Maps/', 'dir')
   mkdir('../3. Final figures/Maps/')
end
if ~exist('../3. Final figures/Spims/', 'dir')
   mkdir('../3. Final figures/Spims/')
end
if ~exist('../3. Final figures/Spectra/', 'dir')
   mkdir('../3. Final figures/Spectra/')
   mkdir('../3. Final figures/Spectra/tex/')
end

addpath('../../algorithms/')
addpath('../../real_data/')

%% Data

disp('Loading data ...')

load('Full2ms.mat'),
load('Full10ms.mat'),

load('../1. Data2use/PCA_NLm.mat'),
load('../1. Data2use/XSNN.mat'),
load('../1. Data2use/XSSS.mat'),
load('../1. Data2use/partial_path.mat')

[n m l] = size(Full2ms);

% Adapt low-quality spectra to fit high-quality ones
Full2ms = Full2ms*5;                   
PCA_NLm = PCA_NLm*5;

% Data structures
Data = {Full2ms, PCA_NLm, Full10ms, XSNN, XSSS};
SampledData = {Full2ms, PCA_NLm, Full10ms, XSNN, XSSS};
ITab = {1:n*m, 1:n*m, 1:n*m, 1:n*m, 1:n*m};

%% output Structures

Nd = length(Data);
B = 5;
FinalSpectra = cell(Nd,B);
FinalMaps = cell(Nd,B);
FinalSpimDisplay = cell(Nd,1);

%% Get normalization at 3 interesting bands
% Full2ms is chosen to prevent saturation.

[ Normalisation ] = getNormalisation( Full2ms );

%% Processing the unmixing

disp('Computing data ...')

% Here, we need to begin the unmixing procedure with an image whose
% abundance maps are clear. This is necessary to easily match the
% next image abundance maps to the reference ones.

ind0 = 2; % reference for matching is PCA_NLm

tab = 1:length(Data);

for i = [ind0 tab(~ismember(tab, ind0))]
    
    if (i==ind0)
        [Spectra,Maps,x0] = perform_unmixing(Data{i}, SampledData{i}, ITab{i});
    else
        [Spectra,Maps,x] = perform_unmixing(Data{i}, SampledData{i}, ITab{i}, x0);
    end
    
    FinalMaps(i,:) = Maps; 
    FinalSpectra(i,:) = transpose(mat2cell(Spectra', ones(1,B),[l])); 
    FinalSpimDisplay{i} = spim2im(Data{i}, Normalisation);
end


%% Output figures

disp('Creating figures ...')

FileDir = '../3. Final figures/Spims/';

% Show spims
cellYlabel = {'Full2ms','PCA_NLm','Full10ms','SNN','SSS'};
Filenames = {'1-X2ms.png','4-X2msDe.png','2-10ms.png','3-SSS.png','6-SNN.png'};

figure,
for ind_y = 1:Nd
    subplot(1,Nd,ind_y)
    imshow(FinalSpimDisplay{ind_y})
    imwrite(FinalSpimDisplay{ind_y},[FileDir Filenames{ind_y}])
    title(cellYlabel{ind_y}),
end

% Partial spim representation
load('../1. Data2use/Partial.mat'), imwrite(spim2im(Partial,Normalisation),[FileDir '5-NoisyPart.png'])

%% Spectra

load('eV.mat')

cellYlabel = {'Full2ms','PCA_NLm','Full10ms','SNN','SSS'};

figure,
for ind_y = 1:Nd
    for ind_x = 1:B
        subplot(B,Nd,ind_y + (ind_x-1)*Nd)
        plot(eV,FinalSpectra{ind_y,ind_x}), grid on,
        if (ind_y == 1) ylabel(['Band #' num2str(ind_x)],'FontSize',18,'FontWeight','Bold'), end
        if (ind_x == 1) title(cellYlabel{ind_y},'FontSize',18,'FontWeight','Bold'), end
    end
end

%% Maps

% Maps Pre-plot
Min = zeros(Nd,B); Max = zeros(Nd,B); 
for ind_x = 1:size(FinalMaps,2), 
    for ind_y = 1:size(FinalMaps,1), 
        Min(ind_y,ind_x) = min(reshape(FinalMaps{ind_y,ind_x},[n*m 1])); 
        Max(ind_y,ind_x) = max(reshape(FinalMaps{ind_y,ind_x},[n*m 1]));
    end; 
end
Min = min(Min,1); Max = max(Max,1);
NormVal = @(x,M,m) (x-m)/(M-m);

NormOne = @(x) (x-min(x(:)))/( max(x(:))-min(x(:)) );

meth = 'AllNormToOne';

cellYlabel = {'Full2ms','PCA_NLm','Full10ms','SNN','SSS'};
figure,
for ind_y = 1:Nd
    for ind_x = 1:B
        
        switch meth
            case 'LetAllAsOutput'
                temp = uint8(255*FinalMaps{ind_y,ind_x});
            case 'AllNormToOne'
                temp = NormOne(FinalMaps{ind_y,ind_x});
            case 'NormPerBand'
                temp = uint8( 255*NormVal(FinalMaps{ind_y,ind_x},Max(ind_x),Min(ind_x) ) );
            otherwise
                temp = uint8(255*FinalMaps{ind_y,ind_x});
        end
        subplot(B,Nd,ind_y + (ind_x-1)*Nd)
        imshow(temp),
        if (ind_y == 1) ylabel(['Band #' num2str(ind_x)],'FontSize',18,'FontWeight','Bold'), end
        if (ind_x == 1) title(cellYlabel{ind_y},'FontSize',18,'FontWeight','Bold'), end
    end
end

% Get chosen bands

% The following input aims at keeping only 3 out of 5 interesting abundance maps.

tab = input('Please give an order tab of bands to keep (e.g. [1 2 3]): ');
FinalMaps2 = FinalMaps(:,tab);
FinalSpectra2 = FinalSpectra(:,tab);
Max = Max(tab); Min = Min(tab);
B2 = length(tab);

% maps
figure,
for ind_y = 1:Nd
    for ind_x = 1:B2
        
        switch meth
            case 'LetAllAsOutput'
                temp = uint8(255*FinalMaps2{ind_y,ind_x});
            case 'AllNormToOne'
                temp = uint8(255*NormOne(FinalMaps2{ind_y,ind_x}));
            case 'NormPerBand'
                temp = uint8( 255*NormVal(FinalMaps2{ind_y,ind_x},Max(ind_x),Min(ind_x) ) );
            otherwise
                temp = uint8(255*FinalMaps2{ind_y,ind_x});
        end
        subplot(B2,Nd,ind_y + (ind_x-1)*Nd)
        imshow(temp),
        imwrite(temp,strcat('../3. Final figures/Maps/',cellYlabel{ind_y},'_',num2str(ind_x),'.png'))
        if (ind_y == 1) ylabel(['Band #' num2str(ind_x)],'FontSize',18,'FontWeight','Bold'), end
        if (ind_x == 1) title(cellYlabel{ind_y},'FontSize',18,'FontWeight','Bold'), end
    end
end

% Plot spectra
Begin = [-500 0 0 ];
End = [3500 5000 4500];
figure,
for ind_x = 1:Nd
    for ind_y = 1:B2
        subplot(B2,Nd,ind_x + (ind_y-1)*Nd)
        plot(eV,FinalSpectra2{ind_x,ind_y}')
        axis([150 900 Begin(ind_y) End(ind_y)])
        grid on
        if (ind_y == 1) title(cellYlabel{ind_x}), end
        if (ind_x == 1) ylabel(['Band #' num2str(ind_y)]), end
    end
end


% save figures

% This part can be commented out if matlab2tikz is installed.
% cf https://github.com/matlab2tikz/matlab2tikz

figure('Position',[10 10 1000 700])
for ind_x = 1:Nd
    for ind_y = 1:B2
        plot(eV,FinalSpectra2{ind_x,ind_y}','Linewidth',4)
        axis([150 900 Begin(ind_y) End(ind_y)])
        set(gca,'FontSize',45)
        grid on
        drawnow
        pause(0.5)
        print(['../3. Final figures/Spectra/spectra_' num2str(ind_y) '_' num2str(ind_x) '.eps' ],'-depsc','-painters')
        %
        % This part can be commented out if matlab2tikz is installed.
        % cf https://github.com/matlab2tikz/matlab2tikz
        % matlab2tikz(['../3. Final figures/Spectra/tex/spectra_' num2str(ind_y) '_' num2str(ind_x) '.tex' ])
    end
end



function [ Normalisation ] = getNormalisation( X )
    % Considering the spectral image X at 3 interesting energy loss, 
    % this function returns for each images a minimum and a scale
    % used for normalization.

    X = transpose(reshape( X,[size(X,1)*size(X,2) size(X,3)] ));

    Bands = [236 346 709];
    R = X(Bands(1),:);
    G = X(Bands(2),:);
    B = X(Bands(3),:);
         

    % third parameter does not exist, so default it to something
    % Normalisation = [ar,br,ag,bg,ab,bb];
    a = @(x) 0.6*min(x(:));
    b = @(x) 1.1*max(x(:))-min(x(:));
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

