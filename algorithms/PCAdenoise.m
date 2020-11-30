function [ Xout ] = PCAdenoise( X, N )
% Denoise a multi-band image by truncature of its PCA at band N
% Input:
%   - X: The 3D multi-band image
%   - N: The order of truncature
% Output:
%   - Xout: The denoised multi-band image

% Get dimensions
[n, m, l] = size(X);

% Reshape and remove mean
Xr = reshape(X,[n*m l]);
Xm = repmat(mean(Xr,1),[n*m 1]);
Xr = Xr-Xm;

% Computes PCA
[V,D] = eig(cov(Xr));
V = fliplr(V);
%keyboard

% Truncature
Vt = V(:,1:N);

% Projection on subspace
Xp = Xr*Vt;

% Output matrix
Xout = reshape( Xp*Vt' + Xm ,[n m l]);



end

