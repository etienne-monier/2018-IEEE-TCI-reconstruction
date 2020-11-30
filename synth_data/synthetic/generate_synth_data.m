clc

% Number of components (Resine, Organic 1, Organic 2 and Calcification)
Nc = 4;

% Abundance maps gaussian blur parameters
kernel_size = 5;
kernel_sig = 3;

% Should the spectra and maps be displayed ?
display_figures = true;

% Sould the variabes be saved ?
save_variables = true;

%
% Abundance maps
%

% All abundance maps were created before. They should be loaded and stored.
Maps = cell(1,Nc);

% Get maps size
map_size = size(rgb2gray(imread('maps/Composante1.png')));

% Resine
Maps{1} = 255*ones(map_size);
% Organic 1
Maps{2} = double(rgb2gray(imread('maps/Composante1.png')));
% Organic 2
Maps{3} = double(rgb2gray(imread('maps/Composante2.png')));
% Calcification
Maps{4} = double(rgb2gray(imread('maps/Composante3.png')));

% Create Gaussian kernel
gauss_kernel = get_blurring_kernel(kernel_size, kernel_sig);

% Proportion for components (Organic1, Organic 2, Calcification).
proportion = [0.5 0.3 0.5];

% Process each map.
for ind = 2:4

    % Get current map
    map = Maps{ind};

    % Thresholding to remove anti-aliasing.
    map(map>0) = 255;

    % Invert white - black.
    map = 255*ones(map_size) - map;

    % Some additional thresholding for Calcification map.
    if (ind==4)
        % center point
        temp = map(30:55, 32:60);
        temp(temp > 250) = 150;
        map(30:55, 32:60) = temp;

        % down left point
        temp = map(70:100, 1:32);
        temp(temp > 250) = 200;
        map(70:100, 1:32) = temp;
    end

    % Perform gaussian blurring.
    map = conv2(map, gauss_kernel, 'same');

    % Set the proportion.
    map = floor(proportion(ind-1)*map);

    % Store result
    Maps{ind} = map;
end

% Create Resina map.
Maps{1} = Maps{1} - Maps{2} - Maps{3} - Maps{4};


% Abundance matrix
P = map_size(1) * map_size(2);
A = zeros([Nc, P]);

% Maps normalization and abundance matrix filling.
for ind=1:Nc
    Maps{ind} = Maps{ind}/255;
    A(ind, :) = Maps{ind}(:);
end


%
% Spectra
%

load('spectra/Spectra.mat');

Spectra = [ResineSpec Orga1Spec Orga2Spec CaSpec];

% Number of chanels
M = length(EvTab);

% Energy loss
eV_offset = 174.4032;
eV_scale = 0.5168;
EnergyLoss = eV_offset + [1:M] * eV_scale;

%
% Linear mixing
%

X_to_be_reshaped = Spectra * A;
X = reshape(X_to_be_reshaped', [map_size, M]);

%
% Save variables to file
%

Elements = {'Resine','Orga1','Orga2','Calcification'};

if save_variables

    DataDescription = 'Synthetic spectrum-image used for IEEEtran TCI paper.';
    
    save '../synthetic_data.mat' X Maps Spectra DataDescription Elements;
end

%
% Display maps and spectra
%

if display_figures

    % Display spectra
    figure
    for ind = 1:4
        subplot(2,2,ind)
        plot(EvTab,Spectra(:, ind))
        title(Elements{ind})
    end

    figure
    for ind = 1:4
        subplot(2,2,ind)
        imshow(uint8(255*Maps{ind}))
        title(Elements{ind})
    end
end


function [kernel] = get_blurring_kernel(k_size, sig)
    % Returns the gaussien blurring kernel based on its size and its
    % standard deviation

    radius = (k_size-1)/2;
    x = - radius : radius;
    [X,Y] = meshgrid(x, x);

    % The gaussian kernel
    kernel = exp(-1 * sqrt(X.^2 + Y.^2) / (2*sig^2));

    % To get sum-to-one kernel.
    normalization = sum(kernel(:));

    % Output kernel
    kernel = kernel/normalization;
end