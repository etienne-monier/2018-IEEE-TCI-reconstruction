function [img_out] = hyper_denoising_nlmeans(img_in,r_search,r_patch,h,sigma)
%% hyper_denoising_nlmeans - performs non-local means denoising
%
% INPUTS
%   img_in   : image to be filtered, should be of type double, of size [m n l]
%   r_search : radio of search window
%   r_patch  : radio of similarity window
%   h        : degree of filtering
%   sigma    : standard deviation of the noise
%
% OUTPUT
%   img_out: the denoised image
%
%  Implementation of the Non local filter proposed for A. Buades, B. Coll and J.M. Morel in
%  "A non-local algorithm for image denoising"
%  Original code by Jose Vicente Manjon Herrera & Antoni Buades
%  Date: 09-03-2006
%
%  Modified to take into account multi-band or hyperspectral 
%  images, and with a different patch similarity measure
%
%  -- By Marie Chabert, Thomas Oberlin & Nicolas Dobigeon, 2017.

%% check inputs
nargin_th = 5;

if nargin<nargin_th || isempty(r_search)
    r_search = 10;
end

if nargin<nargin_th || isempty(r_patch)
    r_patch = 3;
end
   
if nargin<nargin_th || isempty(h)
    h = sigma*0.35;
end

% Size of the image
[m, n, l]=size(img_in);
  
% Memory for the output
img_out=zeros(m,n,l);

% Replicate the boundaries of the input image
[img_in2] = mypadarray(img_in,[r_patch r_patch]);
% img_in2 = padarray(img_in,[r_patch r_patch],'symmetric');
 
bar = waitbar(0,'Please wait...');
for i=1:m 
    for j=1:n 
        % run through all pixels
                 
        i1 = i+ r_patch;
        j1 = j+ r_patch;
        
        % Extract current patch
        W1= img_in2(i1-r_patch:i1+r_patch , j1-r_patch:j1+r_patch,:);
         
        wmax=0; % max similarity
        average=0; % denoised pixel
        sweight=0; % normalization factor
         
        rmin = max(i1-r_search,r_patch+1);
        rmax = min(i1+r_search,m+r_patch);
        smin = max(j1-r_search,r_patch+1);
        smax = min(j1+r_search,n+r_patch);
         
        for r=rmin:1:rmax
            for s=smin:1:smax 
                % run through search window
                                               
                % central pixel: do not compare
                if(r==i1 && s==j1) continue; end;
                
                % extract patch
                W2= img_in2(r-r_patch:r+r_patch , s-r_patch:s+r_patch,:);                
                 
                % computes similarity
                d = sum((W1(:)-W2(:)).^2);
                d = d / l / (2*r_patch+1)^2;
                d = d-2*sigma^2;
                d(d<0) = 0;
                w=exp(-d/h^2);                 
                             
                % updates max, average and normalization factor
                if w>wmax                
                    wmax=w;                   
                end
                sweight = sweight + w;
                average = average + w*img_in2(r,s,:);                                  
            end 
        end
        
        % set w(i,i) again with maximum weight
        average = average + wmax*img_in2(i1,j1,:);
        sweight = sweight + wmax;
                   
        if sweight > 0
            img_out(i,j,:) = average / sweight;
        else
            img_out(i,j,:) = img_in(i,j,:);
        end   
        
        waitbar((j + (i-1)*n) / (n*m))
    end
end

close(bar)
end

function [imout] = mypadarray(im,p_size)

[n,m,l] = size(im);

rx = p_size(1); ry = p_size(2);

imout = zeros(n+2*ry, m+2*rx, l);

% center kept image
imout(ry+[1:n],rx+[1:m],:) = im;
% left and right
imout(ry+[1:n],1:rx,:) = flip(im(:,1:rx,:),2);
imout(ry+[1:n],rx+m+[1:rx],:) = flip(im(:,(m-rx+1):m,:),2);
% up and down
imout(1:ry,rx+[1:m],:) = flip(im(1:ry,:,:),1);
imout(ry+n+[1:ry],rx+[1:m],:) = flip(im((n-ry+1):n,:,:),1);
% up left
imout(1:ry,1:rx,:) = flip(flip( im(1:ry,1:rx,:) ,2),1);
% up right
imout(1:ry,rx+m+[1:rx],:) = flip(flip( im(1:ry,(m-rx+1):m,:) ,2),1);
% down right
imout(ry+n+[1:ry],rx+m+[1:rx],:) = flip(flip( im((n-ry+1):n,(m-rx+1):m,:) ,2),1);
% down left
imout(ry+n+[1:ry],1:rx,:) = flip(flip( im((n-ry+1):n,1:rx,:) ,2),1);


end
