function [ sigma ] = madNoiseEst( X )
% Estimates the noise level of a n*m matrix X. The maximum absolute
%  deviation algorithm is computed for each column of X. The median value of
%  these estimates are then chosen to be an estimation of the noise of X.
%
% Input:
%   X: a n*m matrix
% Output:
%   sigma: the estimated X noise level
%
% Donoho, D.L.; I.M. Johnstone (1994), “Ideal spatial adaptation by wavelet shrinkage,” Biometrika, vol 81, pp. 425–455. 
%
%
% Equivalent wavelet toolbox code
%   [c,l] = wavedec(Xi,2,'db3');
%   sigmai = wnoisest(c,l,1);



% Noise estimation using non-toolbox functions

Nb = size(X,1);
Ns = size(X,2);

% Computes db3 decomposition up to second level.
Jmin = 2;
options.wavelet_type = 'daubechies';
options.wavelet_vm = 3;
J = nextpow2(Nb)-1;
nmin = 2^J;

tmp = zeros(1,Ns);
for ind = 1:Ns
    MW = perform_wavelet_transform(X(1:nmin,ind), J-2, +1, options);
    tmp(ind) = median(abs(MW((2^(J-1)+1):2^(J))))/0.6745;
end
sigma = median(tmp);


end

